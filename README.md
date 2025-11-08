# 🛰️ Drone Operations Optimization — Greedy and Divide-and-Conquer Algorithms

**Contributors:**
- **Algorithm Design & Analysis:** Praneeth Buchepalli  
- **Implementation & Experiments:** Kamal Kandula  

---

## 🧩 Problem Statement

With the rapid adoption of **autonomous drones** for delivery, logistics, and monitoring, modern cities face two major computational challenges:
1. **Efficient Resource Allocation** – Drones must share limited charging infrastructure without overlap or delay.  
2. **Safety & Collision Avoidance** – Drones in dense airspace must be continuously monitored to detect near-miss events in real time.

This project models these real-world challenges and solves them using two fundamental algorithmic paradigms — **Greedy Scheduling** and **Divide-and-Conquer Closest Pair Detection** — ensuring both efficiency and correctness in drone operations.

---

## 📘 About the Project

**Drone Operations Optimization** demonstrates how **classical algorithms** can be applied to modern **urban air mobility** systems.

It introduces two computational modules:

- **Charging Slot Allocation (Greedy Algorithm):**  
  Selects the maximum number of non-overlapping drone charging appointments using an *Earliest-Finish-Time* scheduling rule.

- **Proximity Alert Detection (Divide-and-Conquer Algorithm):**  
  Efficiently identifies the two drones closest to each other among hundreds or thousands of coordinates in real time.

The project combines **algorithmic theory, Python implementation, and experimental validation** to show that classical algorithmic approaches remain powerful for current drone management systems.

---

## ⚙️ Key Features

- **Greedy Charging Scheduler** – Optimally assigns non-overlapping time slots for drone recharging.  
- **Divide-and-Conquer Proximity Detector** – Identifies the closest pair of drones efficiently in O(n log n) time.  
- **Runtime Visualization** – Generates graphs comparing empirical and theoretical complexities.  
- **Experimental Verification** – Benchmarks both algorithms using increasing data sizes.  
- **Result Exports** – Automatically saves runtime data as CSV and figures in `/artifacts/`.  
- **Scalable Design** – Methods are generalizable to multi-hub and multi-drone environments.  

---

## 🧠 Technologies Used

- **Programming Language:** Python 3.13  
- **Libraries:**
  - NumPy – Array and mathematical operations  
  - Matplotlib – Visualization of runtime performance  
  - Pandas – CSV data handling and logging  
- **Development Tools:**
  - Jupyter / VSCode – Development and debugging  
  - GitHub – Version control and code hosting  
- **Operating System:** Windows 11 (tested on Ryzen 7, 16 GB RAM)  

---

## 📚 Algorithms Overview

### 1. **Greedy Scheduling (Charging Slot Allocation)**
- **Problem Type:** Interval Scheduling  
- **Strategy:** Select request with the earliest finishing time  
- **Complexity:** O(n log n) (sorting dominates)  
- **Optimality:** Proven by *exchange argument*  

### 2. **Divide-and-Conquer Closest Pair (Proximity Detection)**
- **Problem Type:** Computational Geometry  
- **Strategy:** Split by median x-coordinate, compute locally, check δ-strip region  
- **Complexity:** O(n log n)  
- **Optimality:** Matches theoretical best-known result by Bentley & Shamos (1976)

---

## 📊 Experimental Results

- **Greedy Scheduler:** Regression slope ≈ 0.998 → matches O(n log n) theory  
- **Divide-and-Conquer:** Regression slope ≈ 0.933 → minor deviation due to cache locality  
- **Runtime Plots:** Saved automatically as `greedy_runtime.png` and `dc_runtime.png`  
- **CSV Outputs:** Detailed timing logs stored in `artifacts/` directory  

## Example outputs:
- Greedy 4000: 0.00075s
- Greedy 16000: 0.00363s
- D&C 1024: 0.00382s
- D&C 4096: 0.01845s

---

## ▶️ How to Run

- pip install numpy matplotlib pandas

- python greedy_dc_experiments.py


