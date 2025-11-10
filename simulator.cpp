#include <bits/stdc++.h>
using namespace std;

//---------------- intilial vals---------------
const int sim_time = 200;        // simulation time in minutes
const int server_nums = 2;          // number of servers
const double avg_service_time = 3;      // mean service time (minutes)
const double avg_interarrival_time = 5; // mean interarrival time (minutes)

//outputs file names
const string STATS_FILE   = "final_summary.csv";
const string TRACE_FILE   = "simulation_trace.txt";
const string EVENTS_FILE  = "event_schedule.csv";
const string SUMMARY_FILE = "final_summary.txt";
//--------------data structures-------------

//record for each customer
struct Customer 
{
    int id;
    int arrival_time;
    int service_start;
    int service_end;
    int server_id;
    int service_time;
    int wait_time;
};
//servers state ie. busy/free etc.
struct Server {
    bool busy = false;
    int totalbusy = 0;
    int current_customer = -1;//bc server is idle (no customer)
    int remaining_time = 0;
    int served = 0;
};

//-----function that generates a rand time interval from an exponential distribution-------
//ensures time is at least 1 minute
int sample_exponential_minutes(double mean, mt19937 &rng) {
    exponential_distribution<double> dist(1.0 / mean);
    return max(1, (int)round(dist(rng)));
}

// -----main loop---------
int main() 
{
    //setup RNG(random num generator)
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    //files initialize here
    ofstream trace(TRACE_FILE);
    ofstream summary(SUMMARY_FILE);
    ofstream events(EVENTS_FILE);
    ofstream stats(STATS_FILE);

    //creates CSV headers
    events << "time,event_type,customer_id,server_id\n";
    stats  << "key,value\n";

    //the simulations var (set up data)
    vector<Server> servers(server_nums);
    vector<Customer> finished;
    queue<Customer> wait_queue;

    int next_arrival = sample_exponential_minutes(avg_interarrival_time, rng);
    int next_customer_id = 1;
    int total_arrivals = 0;

    trace << "Simulation start (" << server_nums << " servers, "
          << sim_time << "min, seed=" <<seed << ")\n";

    // ---------------- MAIN LOOP ----------
    for (int t = 0; t <= sim_time; ++t) {

        //first checking for new arrival
        if (t == next_arrival) {
            Customer c;
            c.arrival_time = t;
            c.id = next_customer_id++;
            c.service_start = -1;
            c.service_end = -1;
            c.wait_time = 0;
            c.server_id = -1;
            c.service_time = sample_exponential_minutes(avg_service_time, rng);
            wait_queue.push(c);
            total_arrivals++;

            events << t << ",arr," << c.id << ",-1\n";
            trace << "[" << t << "] ARRIVAL: customer " << c.id << "\n";

            next_arrival = t + sample_exponential_minutes(avg_interarrival_time, rng);
        }

        //update the busy servers
        for (int i = 0; i < server_nums; ++i) {
            if (servers[i].busy) {
                servers[i].remaining_time--;
                servers[i].totalbusy++;
                if (servers[i].remaining_time <= 0) {
                    servers[i].busy = false;
                    trace << "[" << t << "] DEPARTURE: customer " << servers[i].current_customer
                          << " from server " << i << "\n";
                    events << t << ",dep," << servers[i].current_customer << "," << i << "\n";
                    servers[i].served++;
                    servers[i].current_customer = -1;
                }
            }
        }

//assign the waiting customers in queue to free servers
        for (int i = 0; i < server_nums; ++i) {
            if (!servers[i].busy && !wait_queue.empty()) {
                Customer c = wait_queue.front(); wait_queue.pop();
                c.service_start = t;
                c.wait_time = c.service_start - c.arrival_time;
                c.service_end = c.service_start + c.service_time;
                c.server_id = i;

                servers[i].busy = true;
                servers[i].remaining_time = c.service_time;
                servers[i].current_customer = c.id;

                finished.push_back(c);
                trace << "[" << t << "] START: customer " << c.id
                      << " -> server " << i << " (serv=" << c.service_time
                      << ", wait=" << c.wait_time << ")\n";
            }
        }
    }
//------------ statistics -----------------
    int completed = 0; 
    double total_wait = 0; 
    double total_response = 0;
    double total_service = 0;
    for (auto &c : finished) {
        if (c.service_start >= 0) {
            total_wait += c.wait_time;
            total_service += c.service_time;
            total_response += (c.service_end - c.arrival_time);
            completed++;
        }
    }
    double avg_wait = completed ? total_wait / completed : 0;
    double avg_response = completed ? total_response / completed : 0;
    double avg_service = completed ? total_service / completed : 0;
    double total_util = 0;

    for (auto &s : servers) total_util += (double)s.totalbusy / sim_time;
    double avg_util = total_util / server_nums;

    // -------------resutl--------------------
    stats << "completed," << completed << "\n";
    stats << "arrivals," << total_arrivals << "\n";
    stats << "avg_wait," << avg_wait << "\n";
    stats << "avg_response," << avg_response << "\n";
    stats << "avg_service," << avg_service << "\n";
    stats << "avg_utilization," << avg_util << "\n";

    for (int i = 0; i < server_nums; ++i) {
        stats << "server_" << i << "_util," << (double)servers[i].totalbusy / sim_time << "\n";
        stats << "server_" << i << "_served," << servers[i].served << "\n";
    }

    summary << "Simulation Time: " << sim_time << " minutes\n";
    summary << "Servers: " << server_nums << "\n";
    summary << "Arrivals: " << total_arrivals << "\n";
    summary << "Completed: " << completed << "\n";
    summary << "Average Waiting Time: " << avg_wait << "\n";
    summary << "Average Service Time: " << avg_service << "\n";
    summary << "Average Response Time: " << avg_response << "\n";
    summary << "Average Server Utilization: " << avg_util << "\n\n";
    for (int i = 0; i < server_nums; ++i)
    {
        summary << "Server " << i << ": Util=" << (double)servers[i].totalbusy / sim_time
                << ", Served=" << servers[i].served << "\n";
    }
    trace << "\nSimulation complete.\n";
    trace.close(); summary.close(); events.close(); stats.close();
    cout << "Simulation finished. Files generated:\n"
         << " - " << TRACE_FILE << "\n"
         << " - " << EVENTS_FILE << "\n"
         << " - " << SUMMARY_FILE << "\n"
         << " - " << STATS_FILE << "\n";
}
