#include <iostream>
#include <queue>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <ctime>

using namespace std;

// ============================================================================
// SECTION 1: ENUMS & CONSTANTS
// ============================================================================

enum EventType {
    ARRIVAL,
    DEPARTURE
};

enum TerminationMode {
    BY_SERVED,  // Terminate after maxServed customers
    BY_TIME     // Terminate after horizonT time units
};

// ============================================================================
// SECTION 2: STRUCTS
// ============================================================================

// ----------------------------------------------------------------------------
// EVENT STRUCT
// Owner: ANGGOTA 1
// TODO: Implement operator< for priority_queue (min-heap)
// ----------------------------------------------------------------------------
struct Event {
    EventType type;
    double time;
    int customerID;
    
    bool operator<(const Event& other) const {
    return time > other.time; // reversed for min-heap
    }

};

// ----------------------------------------------------------------------------
// STATE STRUCT
// Owner: ANGGOTA 1 (arrival-related) + ANGGOTA 2 (service-related)
// TODO: Complete state variables
// ----------------------------------------------------------------------------
struct State {
    double clock;                  // Current simulation time
    int numInSystem;               // Total customers in system (queue + service)
    bool serverBusy;               // Server status
    double nextArrivalTime;        // Scheduled next arrival time
    queue<double> arrivalTimes;    // Queue of customer arrival times
    
    State() {
        clock = 0.0;
        numInSystem = 0;
        serverBusy = false;
        nextArrivalTime = 0.0;
        while (!arrivalTimes.empty()) arrivalTimes.pop();
    }

};

// ----------------------------------------------------------------------------
// STATS STRUCT
// Owner: ANGGOTA 2
// TODO: Add tracking variables for statistics
// ----------------------------------------------------------------------------
struct Stats {
    double totalDelay;       // Sum of time in system (W)
    double areaQ;            // Time-weighted queue length integral
    double areaB;            // Time-weighted server busy integral
    int numServed;           // Count of departed customers
    int numArrived;          // Count of arrived customers
    double warmupEndTime;    // Time when warmup period ends // Kapan trigger time buat warmupend nya?
    
    // TODO ANGGOTA 2: Add warmup-related tracking
    Stats() {
        totalDelay = 0.0;
        areaQ = 0.0;
        areaB = 0.0;
        numServed = 0;
        numArrived = 0;
        warmupEndTime = 0.0;
    }
};

// ----------------------------------------------------------------------------
// PARAMS STRUCT
// Owner: ANGGOTA 2
// TODO: Add all simulation parameters
// ----------------------------------------------------------------------------
struct Params {
    double lambda;              // Arrival rate
    double mu;                  // Service rate
    int maxServed;              // Max customers (BY_SERVED mode)
    double horizonT;            // Time horizon (BY_TIME mode)
    int warmup;                 // Warmup period (customers)
    int seed;                   // Random seed
    int queueCap;               // Queue capacity (-1 = unlimited)
    TerminationMode termMode;   // Termination mode
    int reps;                   // Repetitions
    string outdir;              // Output directory
    
    Params() {
        lambda = 0.0;
        mu = 1.0;
        maxServed = 10000;      // yang penting gede
        horizonT = 10000.0;     // yang penting gede
        warmup = 10;            // 10 aja gapapa kayanya?
        seed = static_cast<int>(time(nullptr)); // biar defaultnya beda terus
        queueCap = -1;
        termMode = BY_TIME;
        reps = 10;
        outdir = "./";
    }

    bool queueNotExploding() {
        return (lambda < mu);
    }


};

// ----------------------------------------------------------------------------
// REPLICATION RESULT STRUCT
// Owner: ANGGOTA 2
// TODO: Store results from single replication
// ----------------------------------------------------------------------------
struct RepResult {
    int repID;
    double avgQ;          // Average queue length (L)
    double utilization;   // Server utilization (rho)
    double avgDelay;      // Average time in system (W)
    double avgWait;       // Average time in queue (Wq)
    int numServed;
    double simTime;
    
    RepResult() {
        repID = 0;
        avgQ = 0;
        utilization = 0.0;
        avgDelay = 0.0;
        avgWait = 0.0;
        numServed = 0;
        simTime = 0;
    }
};

// ============================================================================
// SECTION 3: RANDOM NUMBER GENERATOR CLASS
// Owner: ANGGOTA 1
// ============================================================================
class RNG {
private:
    mt19937_64 generator;

public:
    RNG(int seed) {
        generator.seed(seed);
    }

    double exponential(double rate) {
        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator);
        return -log(u) / rate;
    }

    void test(double rate, int samples = 1000) {
        double sum = 0.0;
        for (int i = 0; i < samples; i++)
            sum += exponential(rate);
        cout << "[RNG TEST] Empirical mean: " << (sum / samples) << " | Theoretical: " << (1.0 / rate) << endl;
    }
};

// ============================================================================
// SECTION 4: DES CLASS
// ============================================================================
class DES {
private:
    State state;
    Stats stats;
    Params params;
    RNG rng;
    priority_queue<Event> FEL;  // Future Event List
    int nextCustomerID;
    
public:
    // ------------------------------------------------------------------------
    // CONSTRUCTOR
    // Owner: ALL (call from main)
    // ------------------------------------------------------------------------
    DES(Params p) : params(p), rng(p.seed) {
        nextCustomerID = 1;
    }
    
    // ------------------------------------------------------------------------
    // INIT - Initialize simulation
    // Owner: ANGGOTA 1
    // TODO: Reset state, stats, schedule first arrival
    // ------------------------------------------------------------------------
    void init() {
        state = State(); // reset all state variables
        stats = Stats(); // Anggota 2 will define this struct
        nextCustomerID = 1;

        while (!FEL.empty()) FEL.pop();

        double firstArrival = rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, firstArrival, nextCustomerID);
        state.nextArrivalTime = firstArrival;

        cout << "[INIT] Simulation initialized at t=0.0 (first arrival at t=" << firstArrival << ")" << endl;
}
    
    // ------------------------------------------------------------------------
    // UPDATE TIME INTEGRALS
    // Owner: ANGGOTA 3
    // TODO: Update area under Q(t) and B(t) curves
    // ------------------------------------------------------------------------
    void updateTimeIntegrals(double previousTime) {
        // TODO ANGGOTA 3:
        // 1. Calculate time delta = state.clock - previousTime
        // 2. If past warmup: update stats.areaQ += numInQueue * delta
        // 3. If past warmup: update stats.areaB += (serverBusy ? 1 : 0) * delta
        // Note: numInQueue = numInSystem - (serverBusy ? 1 : 0)
    }
    
    // ------------------------------------------------------------------------
    // HANDLE ARRIVAL
    // Owner: ANGGOTA 1
    // TODO: Process arrival event
    // ------------------------------------------------------------------------
    void handleArrival(Event e) {
        double prevTime = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(prevTime);

        // If server is busy and queue is full, reject
        if (state.serverBusy && params.queueCap != -1 && (int)state.arrivalTimes.size() >= params.queueCap) {
            cout << "[ARRIVAL] Customer " << e.customerID << " rejected (queue full)" << endl;
            return; // stop here, customer rejected
        }

        stats.numArrived++;                                         // kalau direject (queue full) diitung ngga? kalo engga mending pindah bawah kayanya
        state.numInSystem++;

        double nextArrival = state.clock + rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, nextArrival, ++nextCustomerID);
        state.nextArrivalTime = nextArrival;

        // If server idle, start service now
        if (!state.serverBusy) {
            state.serverBusy = true;
            double departTime = state.clock + rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, departTime, e.customerID);
        } else {
            // Server busy → add to queue
            state.arrivalTimes.push(e.time);
            }

            cout << "[ARRIVAL] Customer " << e.customerID <<" arrived at t=" << fixed << setprecision(4) << e.time << endl;
        } 
    }
    
    // ------------------------------------------------------------------------
    // HANDLE DEPARTURE
    // Owner: ANGGOTA 2
    // TODO: Process departure event
    // ------------------------------------------------------------------------
    void handleDeparture(Event e) {
        double prevTime = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(prevTime);

        stats.numServed++;
        
        if (!state.arrivalTimes.empty()) {
            double arrived_at = state.arrivalTimes.front();
            state.arrivalTimes.pop();
            double delay = (e.time - arrived_at);       // artinya delay di kita itungannya dari dateng sampe selesai service?
            if (isWarmupComplete()) {                   // kayanya kata bu novera sampe mulai service deh
                stats.totalDelay += delay;
            }
        }

        state.numInSystem--;

        if(state.numInSystem > 0) {                     // lanjut service antrian selanjutnya
            double serviceTime = rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, state.clock + serviceTime, e.customerID); // customerID nya gimana ya? soalnya kan dia ngga tau customerID siapa yang dilayani selanjutnya
                                                                               // mungkin kita perlu simpan customerID di queue juga?
        } else {
            state.serverBusy = false;                   // servernya nganggur
        }
        
        cout << "[DEPARTURE] Customer " << e.customerID << " at t=" << fixed << setprecision(4) << e.time << endl;
    }
    
    // ------------------------------------------------------------------------
    // CHECK TERMINATION
    // Owner: ANGGOTA 2
    // TODO: Check if simulation should terminate
    // ------------------------------------------------------------------------
    bool shouldTerminate() {
        if (params.termMode == BY_SERVED) {
            return (stats.numServed - params.warmup) >= params.maxServed;   // harusnya yang dateng pas warmup belom diitung
        } else {
            return (state.clock - stats.warmupEndTime)  >= params.horizonT; // harusnya waktu warmup ngga masuk ke horizonT?
        }
        return false; // PLACEHOLDER
    }
    
    // ------------------------------------------------------------------------
    // CHECK WARMUP
    // Owner: ANGGOTA 2
    // TODO: Check if warmup period has ended
    // ------------------------------------------------------------------------
    bool isWarmupComplete() {
        return stats.numServed > params.warmup; 
    }

    bool isWarmupCompletedRightAtTheEndOfThisService() {                    // belom tau sih buat apa
        return stats.numServed == params.warmup;
    }
    
    // ------------------------------------------------------------------------
    // RUN - Main simulation loop
    // Owner: ANGGOTA 1 + ANGGOTA 2 (integration)
    // TODO: Event processing loop
    // ------------------------------------------------------------------------
    RepResult run() {
        // TODO ANGGOTA 1 & 2:
        // 1. Call init()
        // 2. While NOT shouldTerminate():
        //    a. Get next event from FEL (top + pop) // idk ya dari top atau ngga, soalnya kalo liat kodenya di dalem handleArrival(),
                                                     // dia bakal nge push arrival selanjutnya, terus baru ngepush departure buat arrival yang ini
                                                     // regardless kalau dilihat dari time nya duluan orang baru dateng atau yang ini selesai service
                                                     // kayanya harus kita ganti jadi ambil event w/ time terkecil dari FEL
        //    b. Check if warmup complete
        //    c. Process event:
        //       - if ARRIVAL: handleArrival(event)
        //       - if DEPARTURE: handleDeparture(event)
        // 3. Calculate final statistics
        // 4. Return RepResult
        
        init();
        
        // MAIN LOOP HERE
        
        return computeResults();
    }
    
    // ------------------------------------------------------------------------
    // COMPUTE RESULTS
    // Owner: ANGGOTA 3
    // TODO: Calculate performance metrics
    // ------------------------------------------------------------------------
    RepResult computeResults() {
        // TODO ANGGOTA 3:
        // 1. Calculate time period (exclude warmup)
        //    period = clock - warmupEndTime
        // 2. Calculate metrics:
        //    avgQ = areaQ / period
        //    utilization = areaB / period
        //    avgDelay = totalDelay / numServed
        //    avgWait = avgDelay - (1.0 / mu)
        // 3. Create and return RepResult
        
        RepResult res;
        res.repID = 0;
        res.avgQ = 0.0;
        res.utilization = 0.0;
        res.avgDelay = 0.0;
        res.avgWait = 0.0;
        res.numServed = stats.numServed;
        res.simTime = state.clock;
        
        return res;
    }
    
    // ------------------------------------------------------------------------
    // SCHEDULE EVENT (Helper)
    // Owner: ANGGOTA 1
    // ------------------------------------------------------------------------
    void scheduleEvent(EventType type, double time, int customerID) {
        // TODO ANGGOTA 1: Create event and push to FEL
        Event e;
        e.type = type;
        e.time = time;
        e.customerID = customerID;    // Jangan update disini deh, masa setiap call scheduleEvent customerID nya di increment?
                                            // padahal tiap customer ID harus call 2 kali, buat arrival sama departure
                                            // mungkin kita ganti jadi scheduleEvent(type, time, customerID)?
        FEL.push(e);
    }

// ============================================================================
// SECTION 5: MULTI-REPLICATION CONTROLLER
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// RUN MULTIPLE REPLICATIONS
// TODO ANGGOTA 2: Loop over replications
// ----------------------------------------------------------------------------
vector<RepResult> runReplications(Params params, int numReps) {
    
    vector<RepResult> results;
    
    for (int i = 1; i <= numReps; i++) {
        cout << "\n=== REPLICATION " << i << " / " << numReps << " ===" << endl;
        
        Params p = params;
        p.seed = params.seed + i*i;  // Different seed per rep
        
        DES sim(p);
        RepResult res = sim.run();
        res.repID = i;
        results.push_back(res);
        
        cout << "Rep " << i << " complete: AvgQ=" << res.avgQ 
             << " Util=" << res.utilization << endl;
    }
    
    return results;
}

// ============================================================================
// SECTION 6: STATISTICS & CONFIDENCE INTERVAL
// Owner: ANGGOTA 3
// ============================================================================

// ----------------------------------------------------------------------------
// SUMMARY STATISTICS STRUCT
// ----------------------------------------------------------------------------
struct Summary {
    string metric;
    double mean;
    double stdDev;
    double ci_lower;
    double ci_upper;
    double ci_width;
};

// ----------------------------------------------------------------------------
// CALCULATE MEAN
// TODO ANGGOTA 3: Calculate average from vector
// ----------------------------------------------------------------------------
double calculateMean(vector<double>& data) {
    // TODO ANGGOTA 3: Sum all values, divide by count
    double sum = 0.0;
    for (double val : data) {
        sum += val;
    }
    return sum / data.size();
}

// ----------------------------------------------------------------------------
// CALCULATE STANDARD DEVIATION
// TODO ANGGOTA 3: Calculate sample std dev
// ----------------------------------------------------------------------------
double calculateStdDev(vector<double>& data, double mean) {
    // TODO ANGGOTA 3:
    // 1. Calculate sum of squared deviations: Σ(xi - mean)²
    // 2. Divide by (n-1) for sample variance
    // 3. Return sqrt(variance)
    
    double sumSq = 0.0;
    for (double val : data) {
        sumSq += (val - mean) * (val - mean);
    }
    double variance = sumSq / (data.size() - 1);
    return sqrt(variance);
}

// ----------------------------------------------------------------------------
// CALCULATE CONFIDENCE INTERVAL
// TODO ANGGOTA 3: Calculate 95% CI using t-distribution
// ----------------------------------------------------------------------------
Summary calculateCI(vector<double>& data, string metricName) {
    // TODO ANGGOTA 3:
    // 1. Calculate mean
    // 2. Calculate std dev
    // 3. Get t-value for 95% CI (hardcode for common n, or use lookup table)
    //    Example t-values: n=10→2.262, n=15→2.145, n=20→2.093, n=30→2.045
    // 4. Calculate margin of error: t * (stdDev / sqrt(n))
    // 5. CI = [mean - margin, mean + margin]
    
    Summary s;
    s.metric = metricName;
    s.mean = calculateMean(data);
    s.stdDev = calculateStdDev(data, s.mean);
    
    int n = data.size();
    double t_value = 2.262;  // TODO: Use proper t-value based on n
    
    double margin = t_value * (s.stdDev / sqrt(n));
    s.ci_lower = s.mean - margin;
    s.ci_upper = s.mean + margin;
    s.ci_width = 2 * margin;
    
    return s;
}

// ----------------------------------------------------------------------------
// COMPUTE ALL SUMMARIES
// TODO ANGGOTA 3: Generate summary for all metrics
// ----------------------------------------------------------------------------
vector<Summary> computeSummaries(vector<RepResult>& results) {
    // TODO ANGGOTA 3:
    // 1. Extract each metric into separate vector
    //    - avgQ, utilization, avgDelay, avgWait
    // 2. Calculate CI for each metric
    // 3. Return vector of summaries
    
    vector<Summary> summaries;
    
    // Extract avgQ
    vector<double> avgQs;
    for (auto& r : results) {
        avgQs.push_back(r.avgQ);
    }
    summaries.push_back(calculateCI(avgQs, "AvgQ"));
    
    // TODO: Repeat for other metrics (utilization, avgDelay, avgWait)
    
    return summaries;
}

// ============================================================================
// SECTION 7: CSV OUTPUT
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// WRITE PER-REPLICATION CSV
// TODO ANGGOTA 2: Export detailed results
// ----------------------------------------------------------------------------
void writePerRepCSV(vector<RepResult>& results, string filename) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open " << filename << endl;
        return;
    }
    
    // Header
    file << "RepID,AvgQ,Utilization,AvgDelay,AvgWait,NumServed,SimTime\n";
    
    // Data rows
    for (auto& r : results) {
        file << r.repID << ","
             << r.avgQ << ","
             << r.utilization << ","
             << r.avgDelay << ","
             << r.avgWait << ","
             << r.numServed << ","
             << r.simTime << "\n";
    }
    
    file.close();
    cout << "Written: " << filename << endl;
}

// ----------------------------------------------------------------------------
// WRITE SUMMARY CSV
// TODO ANGGOTA 2: Export confidence intervals
// ----------------------------------------------------------------------------
void writeSummaryCSV(vector<Summary>& summaries, string filename) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open " << filename << endl;
        return;
    }
    
    // Header
    file << "Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width\n";
    
    // Data rows
    for (auto& s : summaries) {
        file << s.metric << ","
             << s.mean << ","
             << s.stdDev << ","
             << s.ci_lower << ","
             << s.ci_upper << ","
             << s.ci_width << "\n";
    }
    
    file.close();
    cout << "Written: " << filename << endl;
}

// ============================================================================
// SECTION 8: CLI PARSER
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// PARSE COMMAND LINE ARGUMENTS
// TODO ANGGOTA 2: Extract parameters from argv
// ----------------------------------------------------------------------------
Params parseArguments(int argc, char* argv[]) {
    
    Params p = Params();
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (i+1 < argc) {
            if (arg == "--lambda") {
                p.lambda = stod(argv[++i]);
            }
            else if (arg == "--mu") {
                p.mu = stod(argv[++i]);
            }
            else if (arg == "--horizonT") {
                p.horizonT = stod(argv[++i]);
            }
            else if (arg == "--maxServed") {
                p.maxServed = stoi(argv[++i]);
            }
            else if (arg == "--warmup") {
                p.warmup = stoi(argv[++i]);
            }
            else if (arg == "--reps") {
                p.reps = stoi(argv[++i]);
            }
            else if (arg == "--seed") {
                p.seed = stoi(argv[++i]);
            }
            else if (arg == "--queueCap") {
                p.queueCap = stoi(argv[++i]);
            }
            else if (arg == "--term") {
                if (argv[++i] == "served") {
                    p.termMode = BY_SERVED;
                } else {
                    p.termMode = BY_TIME;
                };
            }
            else if (arg == "--outdir") {
                p.outdir = argv[++i];
            }
        } else if (arg == "--help") {
            printHelp();
            exit(0);
        }
    }
    
    if (!p.queueNotExploding()) {
        cerr << "WARNING: System unstable (lambda >= mu)" << endl;
    }
    
    return p;
}

// ----------------------------------------------------------------------------
// PRINT HELP MESSAGE
// TODO ANGGOTA 2: Display usage information
// ----------------------------------------------------------------------------
void printHelp() {
    // TODO ANGGOTA 2: Print all available parameters and examples
    cout << "Usage: ./des_sim [OPTIONS]\n";
    cout << "\nOptions:\n";
    cout << "  --lambda <value>    Arrival rate (default: 0.0)\n";
    cout << "  --mu <value>        Service rate (default: 1.0)\n";
    cout << "  --horizonT <value>  Time limit before termination (default: 10000.0)\n";
    cout << "  --maxServed <value> N of Served customers before termination (default: 10000)\n";
    cout << "  --warmup <value>    N of starting customers not included in the metrics (default: 10)\n";
    cout << "  --reps <value>      N of repeated simulations (default: 10)\n";
    cout << "  --seed <value>      RNG seed (default: static_cast<int>(time(nullptr)))\n";
    cout << "  --queueCap <value>  Maximum queue length (default: -1 (-1 = no queue limit))\n";
    cout << "  --term <value>      Termination mode (default: time)\n";
    cout << "  --outdir <value>    Output directory (default: ./)\n";
    // TODO: Add all other options
    cout << "\nExample:\n";
    cout << "  ./DES --lambda 0.9 --mu 1.0 --maxServed 20000 --warmup 1000 --reps 10 --seed 12345 --queueCap -1 --term served --outdir ./\n";
}

// ============================================================================
// SECTION 9: VALIDATION & ANALYSIS
// Owner: ANGGOTA 3
// ============================================================================

// ----------------------------------------------------------------------------
// VALIDATE WITH THEORETICAL RESULTS
// TODO ANGGOTA 3: Compare simulation vs M/M/1 theory
// ----------------------------------------------------------------------------
void validateResults(Params p, vector<Summary>& summaries) {
    // TODO ANGGOTA 3:
    // 1. Calculate theoretical M/M/1 results:
    //    rho = lambda / mu
    //    L = rho / (1 - rho)
    //    W = 1 / (mu - lambda)
    //    Wq = rho / (mu - lambda)
    //    Lq = lambda * Wq
    // 2. Compare with simulation results
    // 3. Calculate % deviation
    // 4. Print comparison table
    
    double rho = p.lambda / p.mu;
    double L_theory = rho / (1 - rho);
    double W_theory = 1.0 / (p.mu - p.lambda);
    double Wq_theory = rho / (p.mu - p.lambda);
    
    cout << "\n=== THEORETICAL VS SIMULATION ===" << endl;
    cout << fixed << setprecision(4);
    cout << "Rho (utilization): " << rho << endl;
    cout << "L (theory): " << L_theory << endl;
    // TODO: Print simulation results and deviation
}

// ----------------------------------------------------------------------------
// VERIFY LITTLE'S LAW
// TODO ANGGOTA 3: Check L = lambda * W
// ----------------------------------------------------------------------------
void verifyLittlesLaw(Params p, vector<Summary>& summaries) {
    // TODO ANGGOTA 3:
    // 1. Get L (avgQ) and W (avgDelay) from summaries
    // 2. Calculate L_calculated = lambda * W
    // 3. Compare with L from simulation
    // 4. Print verification result
    
    cout << "\n=== LITTLE'S LAW VERIFICATION ===" << endl;
    // TODO: Implement verification
}

// ============================================================================
// SECTION 10: MAIN FUNCTION
// Owner: ALL (integration point)
// ============================================================================

int main(int argc, char* argv[]) {
    cout << "==================================================" << endl;
    cout << "  M/M/1 Queue Discrete Event Simulation" << endl;
    cout << "==================================================" << endl;
    
    // TODO ANGGOTA 2: Parse command line arguments
    Params params = parseArguments(argc, argv);
    
    // Print configuration
    cout << "\nConfiguration:" << endl;
    cout << "  Lambda: " << params.lambda << endl;
    cout << "  Mu: " << params.mu << endl;
    cout << "  Rho: " << (params.lambda / params.mu) << endl;
    cout << "  Replications: " << params.reps << endl;
    
    // TODO ANGGOTA 2: Run replications
    vector<RepResult> results = runReplications(params, params.reps);
    
    // TODO ANGGOTA 3: Compute summaries
    vector<Summary> summaries = computeSummaries(results);
    
    // TODO ANGGOTA 3: Print results
    cout << "\n=== SUMMARY STATISTICS ===" << endl;
    for (auto& s : summaries) {
        cout << s.metric << ": " << s.mean 
             << " [" << s.ci_lower << ", " << s.ci_upper << "]" << endl;
    }
    
    // TODO ANGGOTA 2: Write CSV files
    writePerRepCSV(results, params.outdir + "results_per_rep.csv");
    writeSummaryCSV(summaries, params.outdir + "summary.csv");
    
    // TODO ANGGOTA 3: Validation
    validateResults(params, summaries);
    verifyLittlesLaw(params, summaries);
    
    cout << "\n=== SIMULATION COMPLETE ===" << endl;
    
    return 0;
}