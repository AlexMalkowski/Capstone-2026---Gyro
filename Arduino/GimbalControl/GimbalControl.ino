/*
 * GimbalControl.ino  --  Scissor-pair gimbal control template
 * Teensy 4.1 + 2x SimpleFOC Mini (DRV8313) + 2x AS5600 (one per I2C bus)
 *
 * Conventions (match the MATLAB sim):
 *   alpha = gimbal A angle, beta = gimbal B angle      [rad]
 *   theta = (alpha + beta)/2      gamma = alpha - beta
 *   alpha = theta + gamma/2       beta  = theta - gamma/2
 *
 * Layers:
 *   1. Axis layer       : per-motor PIDF + feedforward -> Uq (voltage torque mode)
 *   2. Pair layer       : theta/gamma commands -> alpha/beta references
 *   3. Supervisor layer : STUBS for trajectory generation and fall logic
 *
 * Status legend:  [DONE] = ported from the MATLAB sim, not yet validated on
 *                         hardware
 *                 [TODO] = incomplete, needs your input
 *
 * NOT TESTED ON HARDWARE. SimpleFOC calls are written from memory of the
 * library API -- check them against the version you install.
 */

#include <SimpleFOC.h>

// ===================================================================
//  1. CONFIGURATION  (EDIT)
// ===================================================================

// --- 1a. Pins [TODO: fill from the Teensy 4.1 pinout; check SimpleFOC's
//         Teensy page for valid 3-PWM pin groupings] -----------------
#define A_PWM1   -1
#define A_PWM2   -1
#define A_PWM3   -1
#define A_EN     -1
#define B_PWM1   -1
#define B_PWM2   -1
#define B_PWM3   -1
#define B_EN     -1
// I2C: Wire (SDA 18/SCL 19) for A, Wire1 (SDA 17/SCL 16) for B  [verify]

// --- 1b. Motor / driver ---------------------------------------------
const int   POLE_PAIRS     = 0;      // [TODO] not in the datasheet: count magnets / measure
const float V_SUPPLY       = 20.0f;  // [V]
const float UQ_LIMIT       = 6.0f;   // [V] same value as hw.Uq_lim in MATLAB
const float PWM_FREQ_HZ    = 25000;  // [Hz]
const float ALIGN_VOLTAGE  = 2.0f;   // [V] voltage used by initFOC alignment
const float GEAR_RATIO     = 1.0f;   // motor:gimbal (1 = direct drive)
const float AXIS_SIGN_A    = 1.0f;   // flip to -1 if positive motion is backwards
const float AXIS_SIGN_B    = 1.0f;

// --- 1c. Timing: copy the values you MEASURE and re-tune in MATLAB ---
const uint32_t TS_US       = 1000;   // control period [us]
const float    TS          = TS_US*1e-6f;

// --- 1d. Nominal plant (copy from MATLAB printout; used for feedforward)
const float J_NOM    = 6.4e-5f;      // [kg m^2]    [TODO] update after inertia is measured
const float B_NOM    = 0.017f;       // [N m s/rad] B_eff
const float KM_NOM   = 0.056f;       // [N m/V]
float FF_VEL_SCALE   = 1.0f;         // 0..1, matches ctrl.ff_vel
float FF_ACC_SCALE   = 1.0f;         // 0..1, matches ctrl.ff_acc

// --- 1e. Controller gains [TODO: paste from MATLAB "Controller gains"]
const float KP = 30.0f;              // V/rad
const float KI = 100.0f;             // V/(rad s)
const float KD = 0.1f;               // V s/rad
const float TF = 1e-3f;              // derivative filter [s]

// --- 1f. Soft limits [TODO: PLACEHOLDERS -- set from the mechanical stops
//         BEFORE enabling closed loop] -------------------------------
const float DEG2RAD = 0.017453292519943295f;
const float RAD2DEG = 57.29577951308232f;
const float ALPHA_MIN = -90.0f*DEG2RAD, ALPHA_MAX = 90.0f*DEG2RAD;
const float BETA_MIN  = -90.0f*DEG2RAD, BETA_MAX  = 90.0f*DEG2RAD;

// --- 1g. Manual command defaults ------------------------------------
const float THETA_VMAX_DEFAULT = 90.0f*DEG2RAD;   // [rad/s]
const float GAMMA_VMAX_DEFAULT = 90.0f*DEG2RAD;   // [rad/s]

// --- 1h. Logging -----------------------------------------------------
const uint32_t LOG_DECIM = 10;       // print every Nth control tick

// ===================================================================
//  2. TYPES
// ===================================================================

struct Reference { float pos = 0, vel = 0, acc = 0; };   // [rad, rad/s, rad/s^2]

// [DONE] Discrete PID, filtered derivative on measurement, conditional
// integration anti-windup. Same structure as simulate_axis in MATLAB.
struct PIDF {
  float Kp, Ki, Kd, Tf, Ts, ulim;
  float I = 0, D = 0, yPrev = 0;

  void reset(float y0) { I = 0; D = 0; yPrev = y0; }

  float update(float ref, float y, float uff) {
    float e = ref - y;
    D = (Tf*D - Kd*(y - yPrev)) / (Tf + Ts);
    yPrev = y;

    float Itry = I + Ki*Ts*e;
    float uTry = Kp*e + Itry + D + uff;
    float uRaw;
    if (fabsf(uTry) > ulim && ((uTry > 0) == (e > 0))) {
      uRaw = Kp*e + I + D + uff;          // saturated and pushing: hold integrator
    } else {
      I = Itry;  uRaw = uTry;
    }
    return constrain(uRaw, -ulim, ulim);
  }
};

struct Axis {
  const char*      name;
  BLDCMotor*       motor;
  float            sign;
  float            zeroOffset = 0;   // [rad] gimbal-angle zero
  PIDF             pid;
  Reference        ref;              // reference for this gimbal
  float            y = 0;            // measured gimbal angle [rad]
  float            u = 0;            // commanded Uq [V]
  float            uff = 0;          // feedforward part of u [V]
};

enum class SysState  { DISABLED, CLOSED_LOOP, FAULT };
enum class RefSource { MANUAL, TRAJECTORY, FALL_RESPONSE };   // last two are stubs

// ===================================================================
//  3. HARDWARE OBJECTS
// ===================================================================
BLDCMotor          motorA(POLE_PAIRS), motorB(POLE_PAIRS);
BLDCDriver3PWM     driverA(A_PWM1, A_PWM2, A_PWM3, A_EN);
BLDCDriver3PWM     driverB(B_PWM1, B_PWM2, B_PWM3, B_EN);
MagneticSensorI2C  sensorA(AS5600_I2C), sensorB(AS5600_I2C);

Axis axisA, axisB;

// ===================================================================
//  4. STATE
// ===================================================================
SysState  sysState  = SysState::DISABLED;
RefSource refSource = RefSource::MANUAL;

// Manual (serial) command state
struct ManualCmd {
  float thetaTarget = 0, gammaTarget = 0;                  // [rad]
  float thetaVmax = THETA_VMAX_DEFAULT, gammaVmax = GAMMA_VMAX_DEFAULT;
  float thetaRef = 0, gammaRef = 0;                        // slewed references
} manual;

uint32_t tLastUs = 0, tickCount = 0, overruns = 0;
bool     logEnabled = false;
char     faultReason[40] = "";

// ===================================================================
//  5. FUNCTION PROTOTYPES
// ===================================================================

// --- Mode control ---
bool  enableClosedLoop();                       // [DONE] bumpless start from current position
void  disableClosedLoop(const char* reason);    // [DONE] zero output, disable drivers

// --- Manual commands (units: radians, rad/s; serial layer converts from deg) ---
void  commandTheta(float theta);                // [DONE]
void  commandGamma(float gamma);                // [DONE]
void  commandThetaGamma(float theta, float gamma);   // [DONE]
void  setThetaSpeedLimit(float vmax);           // [DONE]
void  setGammaSpeedLimit(float vmax);           // [DONE]
void  holdPosition();                           // [DONE] freeze at current refs
void  zeroAxes();                               // [DONE] current pose becomes alpha = beta = 0

// --- Control tick pieces ---
void  controlTick();                            // [DONE] runs every TS_US
void  readAxes();                               // [DONE]
void  updateReference();                        // [DONE] dispatch on refSource
void  manualReferenceStep(Reference& th, Reference& ga);   // [DONE] velocity-limited slew
void  computeFeedforward(Axis& ax);             // [DONE]
void  checkFaults();                            // [TODO] stub

// --- Serial ---
void  serialPoll();                             // [DONE] non-blocking line reader
void  handleCommand(char* line);                // [DONE] parser (stub entries for traj/fall)
void  printStatus();
void  printHelp();
void  logSample();

// --- Supervisor layer (future larger control model) ---
void  supervisorUpdate();                       // [TODO] mode logic, runs every tick before updateReference
bool  trajectoryStart();                        // [TODO] load / start a trajectory
bool  trajectoryNext(Reference& th, Reference& ga);        // [TODO] one sample of theta/gamma reference
bool  fallDetected();                           // [TODO] detection logic (inputs TBD)
void  fallResponseNext(Reference& th, Reference& ga);      // [TODO] reference during fall response

// ===================================================================
//  6. SETUP / LOOP
// ===================================================================

void setup() {
  Serial.begin(115200);
  uint32_t t0 = millis();
  while (!Serial && millis() - t0 < 2000) {}

  // Refuse to run with unfilled config
  if (POLE_PAIRS <= 0 || A_PWM1 < 0 || B_PWM1 < 0 || A_EN < 0 || B_EN < 0) {
    Serial.println("CONFIG ERROR: set pins and POLE_PAIRS in section 1.");
    while (true) { delay(1000); }
  }

  // Sensors: one AS5600 per bus
  sensorA.init(&Wire);
  sensorB.init(&Wire1);
  // [TODO] confirm the I2C clock (400 kHz) in your SimpleFOC version

  // Drivers
  driverA.voltage_power_supply = V_SUPPLY;  driverA.voltage_limit = UQ_LIMIT;
  driverB.voltage_power_supply = V_SUPPLY;  driverB.voltage_limit = UQ_LIMIT;
  driverA.pwm_frequency = PWM_FREQ_HZ;      driverB.pwm_frequency = PWM_FREQ_HZ;
  driverA.init();  driverB.init();

  // Motors: voltage torque mode, command = Uq [V].
  // Do NOT set phase_resistance: that switches move() to current semantics
  // in newer SimpleFOC versions and would double-count back-EMF.
  motorA.linkSensor(&sensorA);  motorA.linkDriver(&driverA);
  motorB.linkSensor(&sensorB);  motorB.linkDriver(&driverB);
  for (BLDCMotor* m : {&motorA, &motorB}) {
    m->controller        = MotionControlType::torque;
    m->torque_controller = TorqueControlType::voltage;
    m->voltage_limit     = UQ_LIMIT;
    m->voltage_sensor_align = ALIGN_VOLTAGE;
    m->init();
  }
  // initFOC moves the rotor during alignment: keep the flywheel clear.
  // [TODO] store sensor_direction / zero_electric_angle to skip calibration
  motorA.initFOC();
  motorB.initFOC();
  motorA.disable();  motorB.disable();

  // Axes
  axisA.name = "A"; axisA.motor = &motorA; axisA.sign = AXIS_SIGN_A;
  axisB.name = "B"; axisB.motor = &motorB; axisB.sign = AXIS_SIGN_B;
  for (Axis* ax : {&axisA, &axisB}) {
    ax->pid.Kp = KP; ax->pid.Ki = KI; ax->pid.Kd = KD; ax->pid.Tf = TF;
    ax->pid.Ts = TS; ax->pid.ulim = UQ_LIMIT;
  }

  tLastUs = micros();
  Serial.println("Gimbal controller ready (DISABLED). Type 'help'.");
}

void loop() {
  motorA.loopFOC();
  motorB.loopFOC();

  uint32_t now = micros();
  uint32_t lag = now - tLastUs;
  if (lag >= TS_US) {
    if (lag > 2*TS_US) overruns++;
    tLastUs += TS_US;
    if (lag > 5*TS_US) tLastUs = now;      // fell far behind: resync, don't burst
    controlTick();
  }

  motorA.move(axisA.u);                    // zero-order hold between ticks
  motorB.move(axisB.u);

  serialPoll();
}

// ===================================================================
//  7. CONTROL TICK
// ===================================================================

void controlTick() {
  readAxes();
  checkFaults();
  supervisorUpdate();

  if (sysState == SysState::CLOSED_LOOP) {
    updateReference();
    computeFeedforward(axisA);
    computeFeedforward(axisB);
    axisA.u = axisA.pid.update(axisA.ref.pos, axisA.y, axisA.uff);
    axisB.u = axisB.pid.update(axisB.ref.pos, axisB.y, axisB.uff);
  } else {
    axisA.u = 0;  axisB.u = 0;
  }

  tickCount++;
  if (logEnabled && (tickCount % LOG_DECIM == 0)) logSample();
}

void readAxes() {
  // shaft_angle is refreshed inside loopFOC()
  axisA.y = axisA.sign*motorA.shaft_angle/GEAR_RATIO - axisA.zeroOffset;
  axisB.y = axisB.sign*motorB.shaft_angle/GEAR_RATIO - axisB.zeroOffset;
}

void updateReference() {
  Reference th, ga;
  switch (refSource) {
    case RefSource::MANUAL:        manualReferenceStep(th, ga);   break;
    case RefSource::TRAJECTORY:    if (!trajectoryNext(th, ga)) { refSource = RefSource::MANUAL; holdPosition(); manualReferenceStep(th, ga); } break;
    case RefSource::FALL_RESPONSE: fallResponseNext(th, ga);      break;
  }

  // alpha = theta + gamma/2, beta = theta - gamma/2  (linear, so vel/acc map the same way)
  axisA.ref.pos = th.pos + ga.pos/2;  axisA.ref.vel = th.vel + ga.vel/2;  axisA.ref.acc = th.acc + ga.acc/2;
  axisB.ref.pos = th.pos - ga.pos/2;  axisB.ref.vel = th.vel - ga.vel/2;  axisB.ref.acc = th.acc - ga.acc/2;

  // Soft limits on each gimbal [TODO: also zero vel/acc when clamped]
  axisA.ref.pos = constrain(axisA.ref.pos, ALPHA_MIN, ALPHA_MAX);
  axisB.ref.pos = constrain(axisB.ref.pos, BETA_MIN,  BETA_MAX);
}

// Velocity-limited slew toward the manual targets.
// [TODO] add an acceleration limit (the sim's trapezoid) and set ref.acc;
//        until then acceleration feedforward is zero for manual moves.
void manualReferenceStep(Reference& th, Reference& ga) {
  float dTh = constrain(manual.thetaTarget - manual.thetaRef, -manual.thetaVmax*TS, manual.thetaVmax*TS);
  float dGa = constrain(manual.gammaTarget - manual.gammaRef, -manual.gammaVmax*TS, manual.gammaVmax*TS);
  manual.thetaRef += dTh;
  manual.gammaRef += dGa;
  th.pos = manual.thetaRef;  th.vel = dTh/TS;  th.acc = 0;
  ga.pos = manual.gammaRef;  ga.vel = dGa/TS;  ga.acc = 0;
}

// [DONE] u_ff = (J*a*scale + B_eff*v*scale)/Km, from the nominal model
void computeFeedforward(Axis& ax) {
  ax.uff = (FF_ACC_SCALE*J_NOM*ax.ref.acc + FF_VEL_SCALE*B_NOM*ax.ref.vel)/KM_NOM;
}

// [TODO] fault checks that call disableClosedLoop(reason), e.g.:
//   - following error |ref - y| above a limit for N ticks
//   - measured speed above a limit
//   - loop overrun count above a limit
//   - sensor read failure / angle jump
//   - Uq saturated for more than T seconds
void checkFaults() {
}

// ===================================================================
//  8. MODE CONTROL
// ===================================================================

bool enableClosedLoop() {
  if (sysState == SysState::FAULT) { Serial.println("In FAULT: clear with 'dis' first."); return false; }
  readAxes();
  for (Axis* ax : {&axisA, &axisB}) {
    ax->pid.reset(ax->y);
    ax->u = 0;  ax->uff = 0;
  }
  // Bumpless: start from where the gimbals are now
  manual.thetaRef = manual.thetaTarget = (axisA.y + axisB.y)/2;
  manual.gammaRef = manual.gammaTarget =  axisA.y - axisB.y;
  refSource = RefSource::MANUAL;
  tLastUs = micros();
  motorA.enable();  motorB.enable();
  sysState = SysState::CLOSED_LOOP;
  Serial.println("Closed loop ENABLED");
  return true;
}

void disableClosedLoop(const char* reason) {
  sysState = (reason && strncmp(reason, "fault", 5) == 0) ? SysState::FAULT : SysState::DISABLED;
  axisA.u = 0;  axisB.u = 0;
  motorA.move(0);  motorB.move(0);
  motorA.disable();  motorB.disable();
  strncpy(faultReason, reason ? reason : "", sizeof(faultReason) - 1);
  Serial.print("Closed loop DISABLED: ");  Serial.println(reason ? reason : "");
}

// ===================================================================
//  9. MANUAL COMMANDS
// ===================================================================

void commandTheta(float theta) { manual.thetaTarget = theta; refSource = RefSource::MANUAL; }
void commandGamma(float gamma) { manual.gammaTarget = gamma; refSource = RefSource::MANUAL; }
void commandThetaGamma(float theta, float gamma) { commandTheta(theta); commandGamma(gamma); }
void setThetaSpeedLimit(float vmax) { manual.thetaVmax = fabsf(vmax); }
void setGammaSpeedLimit(float vmax) { manual.gammaVmax = fabsf(vmax); }

void holdPosition() {
  manual.thetaTarget = manual.thetaRef;
  manual.gammaTarget = manual.gammaRef;
  refSource = RefSource::MANUAL;
}

void zeroAxes() {
  if (sysState == SysState::CLOSED_LOOP) { Serial.println("Disable closed loop before zeroing."); return; }
  axisA.zeroOffset += axisA.y;
  axisB.zeroOffset += axisB.y;
  readAxes();
  Serial.println("Zeroed: alpha = beta = 0 at current pose");
}

// ===================================================================
//  10. SERIAL INTERFACE   (angles in DEGREES, speeds in DEG/S)
// ===================================================================

void serialPoll() {
  static char buf[64];
  static uint8_t n = 0;
  while (Serial.available()) {
    char c = Serial.read();
    if (c == '\n' || c == '\r') {
      if (n > 0) { buf[n] = 0; handleCommand(buf); n = 0; }
    } else if (n < sizeof(buf) - 1) {
      buf[n++] = tolower(c);
    }
  }
}

void handleCommand(char* line) {
  char* cmd = strtok(line, " ,");
  char* a1  = strtok(NULL, " ,");
  char* a2  = strtok(NULL, " ,");
  if (!cmd) return;

  if      (!strcmp(cmd, "en"))     enableClosedLoop();
  else if (!strcmp(cmd, "dis"))    disableClosedLoop("user request");
  else if (!strcmp(cmd, "th")  && a1)  commandTheta(DEG2RAD*atof(a1));
  else if (!strcmp(cmd, "ga")  && a1)  commandGamma(DEG2RAD*atof(a1));
  else if (!strcmp(cmd, "pos") && a1 && a2) commandThetaGamma(DEG2RAD*atof(a1), DEG2RAD*atof(a2));
  else if (!strcmp(cmd, "vth") && a1)  setThetaSpeedLimit(DEG2RAD*atof(a1));
  else if (!strcmp(cmd, "vga") && a1)  setGammaSpeedLimit(DEG2RAD*atof(a1));
  else if (!strcmp(cmd, "hold"))   holdPosition();
  else if (!strcmp(cmd, "zero"))   zeroAxes();
  else if (!strcmp(cmd, "log"))    logEnabled = (a1 && !strcmp(a1, "on"));
  else if (!strcmp(cmd, "status")) printStatus();
  else if (!strcmp(cmd, "traj"))   { if (!trajectoryStart()) Serial.println("trajectory: not implemented"); }
  else if (!strcmp(cmd, "help"))   printHelp();
  // [TODO] "kp/ki/kd/tf <val>" live gain tuning, "ffv/ffa <0..1>", "fall" test trigger
  else Serial.println("? (type 'help')");
}

void printHelp() {
  Serial.println("en | dis            enable / disable closed loop");
  Serial.println("th <deg>            theta target       ga <deg>   gamma target");
  Serial.println("pos <theta> <gamma> both targets [deg]");
  Serial.println("vth <deg/s>         theta speed limit  vga <deg/s> gamma speed limit");
  Serial.println("hold | zero         freeze refs | set current pose as zero (disabled only)");
  Serial.println("log on|off          CSV stream: t_ms,refA,yA,uA,refB,yB,uB (deg, V)");
  Serial.println("status | traj | help");
}

void printStatus() {
  Serial.print("state=");
  Serial.print(sysState == SysState::CLOSED_LOOP ? "CLOSED_LOOP" : sysState == SysState::FAULT ? "FAULT" : "DISABLED");
  if (sysState == SysState::FAULT) { Serial.print(" ("); Serial.print(faultReason); Serial.print(")"); }
  Serial.print("  alpha="); Serial.print(RAD2DEG*axisA.y, 2);
  Serial.print("  beta=");  Serial.print(RAD2DEG*axisB.y, 2);
  Serial.print("  theta="); Serial.print(RAD2DEG*(axisA.y + axisB.y)/2, 2);
  Serial.print("  gamma="); Serial.print(RAD2DEG*(axisA.y - axisB.y), 2);
  Serial.print("  targets(th,ga)="); Serial.print(RAD2DEG*manual.thetaTarget, 1);
  Serial.print(","); Serial.print(RAD2DEG*manual.gammaTarget, 1);
  Serial.print("  overruns="); Serial.println(overruns);
}

void logSample() {
  Serial.print(millis());                      Serial.print(',');
  Serial.print(RAD2DEG*axisA.ref.pos, 3);      Serial.print(',');
  Serial.print(RAD2DEG*axisA.y, 3);            Serial.print(',');
  Serial.print(axisA.u, 3);                    Serial.print(',');
  Serial.print(RAD2DEG*axisB.ref.pos, 3);      Serial.print(',');
  Serial.print(RAD2DEG*axisB.y, 3);            Serial.print(',');
  Serial.println(axisB.u, 3);
}

// ===================================================================
//  11. SUPERVISOR LAYER  (STUBS for the larger control model)
//  Intended flow, each tick, before updateReference():
//    supervisorUpdate() decides refSource:
//       normal balance/hold  -> MANUAL or TRAJECTORY
//       fallDetected()       -> FALL_RESPONSE (and may override limits)
//  The axis layer never changes: it only sees alpha/beta references.
//  ===================================================================

// [TODO] mode logic, e.g. switch refSource on fallDetected(), handle
// trajectory completion, log mode transitions.
void supervisorUpdate() {
}

// [TODO] start a trajectory. Options: port point_to_point() from the MATLAB
// script and generate online, or stream a precomputed table (export from MATLAB).
bool trajectoryStart() {
  return false;
}

// [TODO] fill th/ga with one sample (pos, vel, acc) per tick; return false when done.
bool trajectoryNext(Reference& th, Reference& ga) {
  (void)th; (void)ga;
  return false;
}

// [TODO] fall detection: inputs and thresholds to be decided (IMU, etc.).
bool fallDetected() {
  return false;
}

// [TODO] reference generation during a fall response.
void fallResponseNext(Reference& th, Reference& ga) {
  (void)th; (void)ga;
}
