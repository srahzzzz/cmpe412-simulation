// Multi-Server Queue Simulation - cleaned and heavily commented
// - Exponential inter-arrival & service times (M/M/c-like behavior)
// - Single-file submission for CMPE412 Term Project 01-02
//
// Outputs:
//   - simulation_trace.txt   : human-readable detailed text trace (minute-by-minute events)
//   - event_schedule.csv     : time_min,event_type,customer_id,server_id
//   - minute_stats.csv       : minute,queue_len,busy_servers,arrivals,departures,avg_wait_so_far,num_in_system,num_in_queue
//   - final_summary.csv      : key,val (final KPIs and per-server lines)
//   - summary.txt            : short human-readable final report (takes info from final_summary.csv)
//
// Usage (compile + run):
//   g++ -std=c++17 -O2 -o simulator simulator_cleaned.cpp
//   ./simulator --servers 2 --simtime 200 --seed 42 --intermean 5.0 --servmean 3.0 --no-finish
//
// Notes about the simulation design :
// - Discrete-minute tick: at each minute we process ALL scheduled events with that minute.
// - Tie-break: DEPARTURE events at minute t are processed before ARRIVAL events at minute t.
// - When a customer starts service we sample an exponential service time and schedule a DEPARTURE event.
// - Interarrival times are sampled from an exponential distribution and converted to an integer minute >= 1.

#include <bits/stdc++.h>
using namespace std;

// -------------------- CONFIGURATION --------------------
static const int DEFAULT_NUM_SERVERS = 2;       // number of servers
static const int DEFAULT_SIM_MINUTES  = 200;   // simulation duration in minutes
static const double DEFAULT_MEAN_INTERARRIVAL = 5.0; // mean interarrival time (minutes)
static const double DEFAULT_MEAN_SERVICE      = 3.0; // mean service time (minutes)

// Output filenames 
const string OUT_TRACE_FILE     = "simulation_trace.txt";   // human-readable event trace
const string OUT_EVENTS_FILE    = "event_schedule.csv";     // event list (csv)
const string OUT_MINUTE_STATS   = "minute_stats.csv";       // per-minute csv stats
const string OUT_SUMMARY_CSV    = "final_summary.csv";      // KPI key,val csv
const string OUT_SUMMARY_TXT    = "summary.txt";            // human-readable summary

// -------------------- TYPES --------------------

// A single customer's canonical record stored for the entire simulation
struct Customer {
    int id = -1;        // unique customer id
    int arrival_min = 0; // minute of arrival
    int service_min = 0;  // sampled service duration in minutes
    int service_start_min = -1; // minute service started (or -1 if not started)
    int service_end_min = -1;   // minute service finished (or -1 if not finished)
    int assigned_server = -1;   // server id assigned (or -1)
};

// A server record
struct Server {
    int id = 0;             // server id (0..num_servers-1)
    bool is_busy = false;   // whether server is currently serving someone
    int service_end_min = -1; // minute when current service ends
    int current_customer = -1; // current customer id or -1
    long long total_busy_minutes = 0; // cumulative busy minutes (utilization numerator)
    int total_served = 0;   // number of customers served by this server
    Server(int i=0): id(i) {}
};

// Event types: arrival or departure
enum EventType { EVENT_DEPARTURE = 0, EVENT_ARRIVAL = 1 };

// An event in the Future Event List (FEL)
struct Event {
    int time_min;      // minute when the event occurs
    EventType type;    // arrival or departure
    int customer_id;   // customer id involved in the event
    int server_id;     // server id for departures; -1 for arrivals
};

// comparator for priority_queue to pop earliest event first
// if times equal, departure events come before arrivals (tie-break rule required by project)
struct EventCompare {
    bool operator()(const Event &a, const Event &b) const {
        if (a.time_min != b.time_min) return a.time_min > b.time_min; // smallest time first
        return a.type > b.type; // EVENT_DEPARTURE (0) before EVENT_ARRIVAL (1)
    }
};

// -------------------- GLOBALS --------------------
int num_servers = DEFAULT_NUM_SERVERS;    // number of parallel servers
int sim_minutes = DEFAULT_SIM_MINUTES;    // total minutes to simulate
unsigned int rng_seed = 0;                // RNG seed
bool rng_seed_provided = false;           // true if seed was passed in CLI
bool finish_pending_after_simtime = true; // if true: finish ongoing services after simtime
double mean_interarrival = DEFAULT_MEAN_INTERARRIVAL; // mean interarrival time
double mean_service = DEFAULT_MEAN_SERVICE;           // mean service time

// Random number generator (mt19937) - seeded in main
static mt19937 rng;

// Future Event List implemented as min-heap priority_queue
priority_queue<Event, vector<Event>, EventCompare> event_queue;

// The canonical collections of servers and customers
vector<Server> server_list;
unordered_map<int, Customer> customers; // key: customer id -> Customer
queue<int> waiting_queue; // stores only customer ids (keeps canonical record authoritative)

// -------------------- STATISTICS (aggregators) --------------------
long long total_wait_minutes = 0;   // sum of waiting times for completed customers (Wq numerator)
long long total_service_minutes = 0; // sum of service times for completed customers
long long total_response_minutes = 0; // sum of response times (end - arrival) for completed customers
long long completed_customers = 0;  // number of customers who completed service
long long started_customers = 0;    // number of customers that started service
long long arrival_count = 0;        // number of arrivals that occurred

int max_queue_length = 0;          // maximum queue length observed
long long customers_who_waited = 0; // number of customers who experienced waiting (start > arrival)

// sums for time-average computations (L and Lq)
long long sum_queue_lengths = 0;  // sum over ticks of queue length (for Lq average)
long long sum_in_system = 0;      // sum over ticks of number in system (L)
int ticks_recorded = 0;           // number of recorded ticks (should be sim_minutes if collecting per tick)

int next_customer_id = 1;         // incremental customer id generator

// file streams for outputs
ofstream out_trace;   // simulation_trace.txt
ofstream out_events;  // event_schedule.csv
ofstream out_minute;  // minute_stats.csv
ofstream out_summary; // final_summary.csv
ofstream out_summary_txt; // summary.txt (human readable)

// -------------------- UTILITIES (documented) --------------------

// Sample an exponential distribution with mean "mean" and return an integer number of minutes >= 1.
// We sample continuous double from exponential(lambda=1/mean), then map to integer minutes by floor.
static int sample_exponential_minutes(double mean) {
    if (!(mean > 0.0)) return 1; // fallback safety
    exponential_distribution<double> dist(1.0 / mean); // lambda = 1/mean
    double sample = dist(rng); // sample in continuous minutes
    int minutes = max(1, (int)floor(sample + 1e-9)); // integer minute at least 1
    return minutes;
}

// Open all output files and write csv headers. Clear comments explain every file header.
static void open_output_files() {
    out_trace.open(OUT_TRACE_FILE);
    out_events.open(OUT_EVENTS_FILE);
    out_minute.open(OUT_MINUTE_STATS);
    out_summary.open(OUT_SUMMARY_CSV);
    out_summary_txt.open(OUT_SUMMARY_TXT);

    // Write headers for event_schedule.csv
    if (out_events.is_open()) out_events << "time_min,event_type,customer_id,server_id\n";

    // Write headers for minute_stats.csv
    if (out_minute.is_open()) out_minute << "minute,queue_len,busy_servers,arrivals,departures,avg_wait_so_far,num_in_system,num_in_queue\n";

    // final_summary.csv header (will append key,value pairs later)
    if (out_summary.is_open()) out_summary << "key,val\n";

    // summary.txt: write top header; human readable body appended at the end
    if (out_summary_txt.is_open()) {
        out_summary_txt << "SIMULATION SUMMARY\n";
        out_summary_txt << "------------------\n";
    }
}

// Close files safely
static void close_output_files() {
    if (out_trace.is_open()) out_trace.close();
    if (out_events.is_open()) out_events.close();
    if (out_minute.is_open()) out_minute.close();
    if (out_summary.is_open()) out_summary.close();
    if (out_summary_txt.is_open()) out_summary_txt.close();
}

// Log a single line both to stdout and to the trace file (trace file is human readable timeline)
static void trace_log(const string &message, int minute) {
    string line = "[" + to_string(minute) + "] " + message;
    cout << line << "\n";
    if (out_trace.is_open()) out_trace << line << "\n";
}

// Write a single event row to event_schedule.csv. event_type is "arr" or "dep".
static void write_event_row(const Event &ev) {
    if (!out_events.is_open()) return;
    string type_str = (ev.type == EVENT_ARRIVAL) ? "arr" : "dep";
    out_events << ev.time_min << "," << type_str << "," << ev.customer_id << "," << ev.server_id << "\n";
}

// Find an index of any free server; return -1 if none available.
static int find_free_server() {
    for (int i = 0; i < (int)server_list.size(); ++i) {
        if (!server_list[i].is_busy) return i;
    }
    return -1;
}

// Assign a given customer (by id) to an available server at 'current_minute'.
// If no server is free, enqueue the customer's id in waiting_queue.
static void assign_customer_to_server(int customer_id, int current_minute) {
    int server_id = find_free_server();
    if (server_id == -1) {
        // no free server: put customer in waiting queue
        waiting_queue.push(customer_id);
        if ((int)waiting_queue.size() > max_queue_length) max_queue_length = (int)waiting_queue.size();
        return;
    }

    // Fetch canonical records
    Customer &c = customers[customer_id];
    Server &s = server_list[server_id];

    // Record service start minute and assigned server
    c.service_start_min = current_minute;
    c.assigned_server = server_id;

    // Sample service duration now (service time is determined at start)
    int service_duration = sample_exponential_minutes(mean_service);
    c.service_min = service_duration;
    c.service_end_min = c.service_start_min + c.service_min;

    // Update server status and busy-time accumulator
    s.is_busy = true;
    s.current_customer = customer_id;
    s.service_end_min = c.service_end_min;
    s.total_busy_minutes += c.service_min; // accumulate busy minutes for utilization

    // Schedule departure event in FEL
    Event dep_ev; dep_ev.time_min = c.service_end_min; dep_ev.type = EVENT_DEPARTURE; dep_ev.customer_id = customer_id; dep_ev.server_id = server_id;
    event_queue.push(dep_ev);
    write_event_row(dep_ev);

    // Stats bookkeeping
    started_customers++;
    if (c.service_start_min > c.arrival_min) customers_who_waited++;

    // Trace log line
    trace_log("ASSIGN: customer " + to_string(customer_id) + " -> server " + to_string(server_id) +
              " (service=" + to_string(c.service_min) + ", end=" + to_string(c.service_end_min) + ")",
              current_minute);
}

// While there is a free server and customers are waiting, assign them in FIFO order
static void drain_waiting_queue(int current_minute) {
    while (!waiting_queue.empty()) {
        int server_id = find_free_server();
        if (server_id == -1) break; // no free servers left
        int customer_id = waiting_queue.front(); waiting_queue.pop();
        assign_customer_to_server(customer_id, current_minute);
    }
}

// Handle departure event: free server, compute metrics for finished customer
static void process_departure(const Event &ev, int current_minute) {
    int server_id = ev.server_id;
    int customer_id = ev.customer_id;

    // defensive server update
    if (server_id >= 0 && server_id < (int)server_list.size()) {
        Server &s = server_list[server_id];
        s.is_busy = false;
        s.current_customer = -1;
        s.service_end_min = -1;
        s.total_served += 1;
    }

    // update customer metrics and global aggregators
    auto it = customers.find(customer_id);
    if (it != customers.end()) {
        Customer &c = it->second;
        c.service_end_min = ev.time_min; // should equal scheduled end

        long long wait_minutes = (long long)c.service_start_min - c.arrival_min; // waiting time
        long long service_minutes = (long long)c.service_min;
        long long response_minutes = (long long)c.service_end_min - c.arrival_min;

        total_wait_minutes += wait_minutes;
        total_service_minutes += service_minutes;
        total_response_minutes += response_minutes;
        completed_customers += 1;

        trace_log("DEPART: customer " + to_string(customer_id) + " from server " + to_string(server_id)
                  + " (wait=" + to_string(wait_minutes) + ", service=" + to_string(service_minutes) + ")",
                  current_minute);
    } else {
        // shouldn't happen, but log for debugging
        trace_log("DEPART: customer not found id=" + to_string(customer_id), current_minute);
    }
}

// Handle arrival event: create canonical customer record (no assignment here)
static void process_arrival(const Event &ev, int current_minute) {
    Customer c;
    c.id = ev.customer_id;
    c.arrival_min = ev.time_min;
    c.service_start_min = -1; c.service_end_min = -1; c.service_min = 0; c.assigned_server = -1;

    customers[c.id] = c;
    arrival_count++;

    trace_log("ARRIVE: customer " + to_string(ev.customer_id), current_minute);
}

// -------------------- MAIN SIMULATION LOOP --------------------
int main(int argc, char **argv) {
    // CLI parsing (basic; unchanged behavior)
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--servers" && i+1 < argc) { num_servers = stoi(argv[++i]); continue; }
        if (arg == "--simtime" && i+1 < argc) { sim_minutes = stoi(argv[++i]); continue; }
        if (arg == "--seed" && i+1 < argc) { rng_seed = (unsigned)stoul(argv[++i]); rng_seed_provided = true; continue; }
        if (arg == "--intermean" && i+1 < argc) { mean_interarrival = stod(argv[++i]); continue; }
        if (arg == "--servmean" && i+1 < argc) { mean_service = stod(argv[++i]); continue; }
        if (arg == "--no-finish") { finish_pending_after_simtime = false; continue; }
        if (arg == "--help") {
            cout << "Usage: " << argv[0] << " [--servers N] [--simtime M] [--seed S] [--intermean x] [--servmean y] [--no-finish]\n";
            return 0;
        }
    }

    // validate CLI values and apply defaults if needed
    if (num_servers <= 0) num_servers = DEFAULT_NUM_SERVERS;
    if (sim_minutes <= 0) sim_minutes = DEFAULT_SIM_MINUTES;

    // seed RNG
    if (!rng_seed_provided) rng_seed = (unsigned)chrono::high_resolution_clock::now().time_since_epoch().count();
    rng.seed(rng_seed);

    // open outputs and write header lines
    open_output_files();

    // trace start line
    trace_log("SIM START: servers=" + to_string(num_servers)
              + " sim_minutes=" + to_string(sim_minutes)
              + " interMean=" + to_string(mean_interarrival)
              + " servMean=" + to_string(mean_service)
              + " seed=" + to_string(rng_seed), 0);

    // initialize server list
    server_list.clear();
    for (int i = 0; i < num_servers; ++i) server_list.emplace_back(i);

    // schedule first arrival
    int first_delta = sample_exponential_minutes(mean_interarrival);
    int first_time = first_delta; // first arrival minute
    Event first_event; first_event.time_min = first_time; first_event.type = EVENT_ARRIVAL; first_event.customer_id = next_customer_id++; first_event.server_id = -1;
    event_queue.push(first_event);
    write_event_row(first_event);

    // main discrete-minute loop
    int minute = 0;
    int last_processed_minute = 0;
    const int SAFETY_MAX_MINUTES = 10'000'000; // safety cap in case of logic errors

    while (true) {
        bool simtime_exceeded = (minute >= sim_minutes);

        // detect whether any server is busy
        bool any_server_busy = false;
        for (auto &s : server_list) if (s.is_busy) { any_server_busy = true; break; }

        bool fel_not_empty = !event_queue.empty();
        bool waiting_not_empty = !waiting_queue.empty();

        // stop condition
        if (simtime_exceeded) {
            if (!finish_pending_after_simtime) {
                // immediate stop at sim_minutes
                break;
            } else {
                // finish pending services (stop only when FEL empty, waiting queue empty, and all servers idle)
                if (!fel_not_empty && !waiting_not_empty && !any_server_busy) break;
            }
        }

        int arrivals_this_min = 0;
        int departures_this_min = 0;

        // 1) Process all departures scheduled at this minute (FEL pops DEPs first due to comparator)
        while (!event_queue.empty()) {
            Event ev = event_queue.top();
            if (ev.time_min != minute) break;
            if (ev.type == EVENT_DEPARTURE) {
                event_queue.pop();
                process_departure(ev, minute);
                departures_this_min++;
                write_event_row(ev);
            } else break; // stop at first arrival - we handle arrivals after departures
        }

        // 2) After departures, assign waiting customers (freed servers now available)
        drain_waiting_queue(minute);

        // 3) Collect arrival events scheduled at this minute (we first collect them so we can schedule next arrivals deterministically)
        vector<Event> arrivals_now;
        while (!event_queue.empty()) {
            Event ev = event_queue.top();
            if (ev.time_min != minute) break;
            if (ev.type == EVENT_ARRIVAL) {
                event_queue.pop();
                arrivals_now.push_back(ev);
            } else break;
        }

        // process each arrival and schedule next arrival
        for (auto &aev : arrivals_now) {
            process_arrival(aev, minute);
            arrivals_this_min++;

            // schedule next arrival based on interarrival distribution
            int next_delta = sample_exponential_minutes(mean_interarrival);
            int next_time = aev.time_min + next_delta;

            if (!simtime_exceeded) {
                // during sim time we schedule next arrivals
                if (next_time <= sim_minutes || finish_pending_after_simtime) {
                    Event nev; nev.time_min = next_time; nev.type = EVENT_ARRIVAL; nev.customer_id = next_customer_id++; nev.server_id = -1;
                    event_queue.push(nev);
                    write_event_row(nev);
                }
            }
        }

        // 4) After creating arrivals, assign newly arrived customers (in increasing customer id to preserve arrival order)
        vector<int> newly_arrived_ids;
        newly_arrived_ids.reserve(arrivals_now.size());
        for (auto &p : customers) {
            const Customer &c = p.second;
            if (c.arrival_min == minute && c.service_start_min == -1) newly_arrived_ids.push_back(c.id);
        }
        sort(newly_arrived_ids.begin(), newly_arrived_ids.end());
        for (int cid : newly_arrived_ids) assign_customer_to_server(cid, minute);

        // try draining waiting queue again (in case some assignments freed servers)
        drain_waiting_queue(minute);

        // 5) Collect per-minute statistics and append to minute_stats.csv
        int queue_len = (int)waiting_queue.size();
        int busy_servers = 0;
        for (auto &s : server_list) if (s.is_busy) busy_servers++;

        int in_system = queue_len + busy_servers; // number of customers present in the system this minute
        sum_queue_lengths += queue_len;
        sum_in_system += in_system;
        ticks_recorded++;

        if (queue_len > max_queue_length) max_queue_length = queue_len;

        double avg_wait_so_far = (completed_customers > 0) ? ((double)total_wait_minutes / (double)completed_customers) : 0.0;
        if (out_minute.is_open()) {
            out_minute << minute << "," << queue_len << "," << busy_servers << "," << arrivals_this_min << "," << departures_this_min
                       << "," << fixed << setprecision(6) << avg_wait_so_far << "," << in_system << "," << queue_len << "\n";
        }

        // compact trace line for readability
        trace_log("TICK: min=" + to_string(minute)
                  + " q=" + to_string(queue_len)
                  + " busy=" + to_string(busy_servers)
                  + " arr=" + to_string(arrivals_this_min)
                  + " dep=" + to_string(departures_this_min)
                  + " avgW=" + to_string(avg_wait_so_far), minute);

        minute++;
        last_processed_minute = minute;

        if (minute > SAFETY_MAX_MINUTES) {
            trace_log("SAFETY STOP: minute exceeded safety cap", minute);
            break;
        }
    } // end main loop

    // -------------------- COMPUTE FINAL KPIS --------------------
    trace_log("SIM END: computing KPIs", last_processed_minute);

    double avg_wait = (completed_customers > 0) ? (double)total_wait_minutes / (double)completed_customers : 0.0; // Wq
    double avg_service = (completed_customers > 0) ? (double)total_service_minutes / (double)completed_customers : 0.0; // S
    double avg_response = (completed_customers > 0) ? (double)total_response_minutes / (double)completed_customers : 0.0; // W
    double throughput = (last_processed_minute > 0) ? (double)completed_customers / (double)last_processed_minute : 0.0; // cust per minute
    double Lq_avg = (ticks_recorded > 0) ? (double)sum_queue_lengths / (double)ticks_recorded : 0.0;
    double L_avg  = (ticks_recorded > 0) ? (double)sum_in_system / (double)ticks_recorded : 0.0;
    double pct_waited = (arrival_count > 0) ? (double)customers_who_waited / (double)arrival_count * 100.0 : 0.0;

    // per-server utilization
    vector<double> server_util;
    double avg_util = 0.0;
    for (auto &s : server_list) {
        double u = (last_processed_minute > 0) ? (double)s.total_busy_minutes / (double)last_processed_minute : 0.0;
        server_util.push_back(u);
        avg_util += u;
    }
    if (!server_util.empty()) avg_util /= (double)server_util.size();

    // write final_summary.csv (key,val)
    if (out_summary.is_open()) {
        out_summary << "completed," << completed_customers << "\n";
        out_summary << "arrivals," << arrival_count << "\n";
        out_summary << "Wq," << avg_wait << "\n";
        out_summary << "S," << avg_service << "\n";
        out_summary << "W," << avg_response << "\n";
        out_summary << "throughput_per_min," << throughput << "\n";
        out_summary << "Lq_avg," << Lq_avg << "\n";
        out_summary << "L_avg," << L_avg << "\n";
        out_summary << "max_q," << max_queue_length << "\n";
        out_summary << "pct_waited," << pct_waited << "\n";
        out_summary << "util_avg," << avg_util << "\n";
        for (int i = 0; i < (int)server_list.size(); ++i) {
            out_summary << "server" << i << "_util," << server_util[i] << "\n";
            out_summary << "server" << i << "_busy," << server_list[i].total_busy_minutes << "\n";
            out_summary << "server" << i << "_served," << server_list[i].total_served << "\n";
        }
    }

    // write summary.txt (human-readable nicely formatted)
    if (out_summary_txt.is_open()) {
        out_summary_txt << "Simulation time (minutes): " << last_processed_minute << "\n";
        out_summary_txt << "Number of servers: " << num_servers << "\n";
        out_summary_txt << "Total arrivals: " << arrival_count << "\n";
        out_summary_txt << "Total completed customers: " << completed_customers << "\n\n";
        out_summary_txt << "Key performance indicators:\n";
        out_summary_txt << " - Average waiting time (Wq): " << avg_wait << " minutes\n";
        out_summary_txt << " - Average service time (S): " << avg_service << " minutes\n";
        out_summary_txt << " - Average response time (W): " << avg_response << " minutes\n";
        out_summary_txt << " - Throughput (cust/min): " << throughput << "\n";
        out_summary_txt << " - Average number in queue (Lq): " << Lq_avg << "\n";
        out_summary_txt << " - Average number in system (L): " << L_avg << "\n";
        out_summary_txt << " - Maximum observed queue length: " << max_queue_length << "\n";
        out_summary_txt << " - Percent customers who waited: " << pct_waited << "%\n";
        out_summary_txt << " - Average server utilization: " << avg_util << "\n\n";
        out_summary_txt << "Per-server details:\n";
        for (int i = 0; i < (int)server_list.size(); ++i) {
            out_summary_txt << " Server " << i << ": util=" << server_util[i]
                            << "  busy_minutes=" << server_list[i].total_busy_minutes
                            << "  served=" << server_list[i].total_served << "\n";
        }
        out_summary_txt << "\nNotes:\n";
        out_summary_txt << " - Use --seed N to reproduce runs (same seed => same outputs)\n";
        out_summary_txt << " - To perform experiments for report, run multiple sims with varying --intermean, --servmean, and --servers\n";
    }

    // final trace output (human-readable)
    trace_log("=== FINAL SUMMARY ===", last_processed_minute);
    trace_log("completed customers: " + to_string(completed_customers), last_processed_minute);
    trace_log("arrivals total: " + to_string(arrival_count), last_processed_minute);
    trace_log("Wq (avg wait): " + to_string(avg_wait), last_processed_minute);
    trace_log("S (avg serv): " + to_string(avg_service), last_processed_minute);
    trace_log("W (avg resp): " + to_string(avg_response), last_processed_minute);
    trace_log("throughput (cust/min): " + to_string(throughput), last_processed_minute);
    trace_log("Lq_avg: " + to_string(Lq_avg) + "   L_avg: " + to_string(L_avg), last_processed_minute);
    trace_log("max queue length: " + to_string(max_queue_length), last_processed_minute);
    trace_log("percent customers who waited: " + to_string(pct_waited) + "%", last_processed_minute);
    trace_log("avg server util: " + to_string(avg_util), last_processed_minute);

    for (int i = 0; i < (int)server_list.size(); ++i) {
        trace_log("server " + to_string(i) + " busy_minutes=" + to_string(server_list[i].total_busy_minutes)
                  + " served=" + to_string(server_list[i].total_served)
                  + " util=" + to_string(server_util[i]), last_processed_minute);
    }

    close_output_files();

    // notify generated files
    cout << "\nGenerated files:\n"
         << " - " << OUT_TRACE_FILE << "\n"
         << " - " << OUT_EVENTS_FILE << "\n"
         << " - " << OUT_MINUTE_STATS << "\n"
         << " - " << OUT_SUMMARY_CSV << "\n"
         << " - " << OUT_SUMMARY_TXT << "\n\n";

    cout << "Tips:\n"
         << "- Use --seed N to reproduce runs (same seed => same outputs).\n"
         << "- Run multiple times with different seeds or parameters for experiments.\n";

    return 0;
}