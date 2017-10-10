# **Unscented Kalman Filter Project: Sensor Fusion**
#### _Lorenzo's version_
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)


[//]: # (Image References)
[results]: ./results.png

---

The goals / steps of this project are the following:
* Utilize an **unscented** kalman filter to estimate the state of a moving object of interest with noisy **laser** and **radar** measurements
* leverage the unscented kalman filter in order to model a **CTRV** (Constant Turn Rate and Velocity magnitude model) process, which state is composed of:
  - x coordinate: *p<sub>x</sub>*
  - y coordinate: *p<sub>y</sub>*
  - velocity: *v*
  - turn-rate: *ψ*
  - change in turn-rate: *ψ<sub>dot</sub>*
* Use sensor fusion to incorporate the information from both lidar and radar sensors.
* Measure the efficiency of the Kalman Filter using Root Mean Square Error (RMSE)

The essential part of the code for this project is the following files:
- **src/ukf.cpp**
- **src/tools.cpp**

After implementing the above steps and passing the algorithm through the [Udacity Simulator](https://github.com/udacity/self-driving-car-sim/releases), we get the below results:

![alt text][results]


---

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./UnscentededKF `

## Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
