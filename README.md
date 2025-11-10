# CMPE412-Simulation Project

A simple single-queue multi-server discrete-event simulator written in C++.

This repository contains a compact simulation that models customer arrivals and service on multiple identical servers. Arrivals and service times are sampled from exponential distributions. The simulator runs in discrete 1-minute steps and produces a trace, an event schedule (CSV), and summary statistics.

## Repository contents

- `simulator.cpp` — The complete simulator source code (single-file). Modify the constants at the top to configure the run.
- `simulator` — Compiled executable (binary produced when you build `simulator.cpp`).
- `event_schedule.csv` — CSV timeline of arrival (arr) and departure (dep) events produced by the last run.
- `simulation_trace.txt` — Human-readable trace of arrivals, service starts, and departures.
- `final_summary.txt` — Plain-text summary with key metrics.
- `final_summary.csv` — CSV of summary statistics (two columns: key,value).

## Overview of the simulator behavior

- Time is discrete and counted in minutes from `0` to `sim_time` (inclusive).
- Customer interarrival and service times are sampled from exponential distributions (rounded, minimum 1 minute).
- Arriving customers enter a FIFO wait queue and are assigned to any free server in order.
- Each server maintains simple counters for total busy time and customers served.
- At the end of the run the simulator writes:
  - `simulation_trace.txt` — detailed timeline useful for debugging.
  - `event_schedule.csv` — easy-to-parse event list for analysis (columns: time,event_type,customer_id,server_id).
  - `final_summary.txt` and `final_summary.csv` — derived statistics (completed arrivals, average wait, response, service time, server utilizations).

## Configuration (edit `simulator.cpp` constants)

At the top of `simulator.cpp` you'll find constants you can change before building:

- `sim_time` — total simulation length in minutes (default 200).
- `server_nums` — number of identical servers (default 2).
- `avg_service_time` — mean of the exponential service time (minutes, default 3).
- `avg_interarrival_time` — mean of the exponential interarrival time (minutes, default 5).

These are simple compile-time constants in the sample; change them and rebuild to run different scenarios.

## Build

The simulator is a single C++ file and can be built with g++ (C++17 or later recommended). Example:

```bash
# from repository root
g++ -std=c++17 -O2 -o simulator simulator.cpp
```
The workspace also contains a VS Code task (`C/C++: g++ build active file`) that will compile the active file in the editor.


## Extending the simulator

This code is intentionally compact. Ideas for extension:

- Add command-line arguments to override constants (e.g., `--sim-time`, `--servers`, `--seed`).
- Replace discrete-time stepping with a full event-driven engine (scheduling the next event directly rather than stepping each minute).
- Add distributions other than exponential, or allow passing a seed for reproducible experiments.
- Collect time-series metrics (e.g., queue length vs time) and output them for plotting.
- If you want reproducible runs, modify the seed initialization in `simulator.cpp` to use a fixed seed instead of time-based seeding.

## Development notes

- The simulator uses the C++ `<random>` library with `exponential_distribution` and `mt19937` RNG.
- Minimum sampled time is clamped to 1 minute (`max(1, round(...))`) to avoid zero-length intervals.
- Outputs are plain text and CSV for ease of parsing in Python, R, or Excel.

## Created by Sarah Nauman and Afagh Izadi Dakhrabadi.
