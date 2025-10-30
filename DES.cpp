#include <iostream>
#include <queue>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

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
    double warmupEndTime;    // Time when warmup period ends
    
    // TODO ANGGOTA 2: Add warmup-related tracking
    // TODO ANGGOTA 2: Initialize in constructor
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
    string outdir;              // Output directory
    
    // TODO ANGGOTA 2: Add constructor with default values
    // TODO ANGGOTA 2: Add validation method (lambda < mu, etc)
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
    
    // TODO ANGGOTA 2: Add constructor
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
        cout << "[RNG TEST] Empirical mean: " << (sum / samples)
             << " | Theoretical: " << (1.0 / rate) << endl;
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
        scheduleEvent(ARRIVAL, firstArrival);
        state.nextArrivalTime = firstArrival;

        cout << "[INIT] Simulation initialized at t=0.0 (first arrival at t=" 
             << firstArrival << ")" << endl;
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

        stats.numArrived++;

        double nextArrival = state.clock + rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, nextArrival);
        state.nextArrivalTime = nextArrival;

        if (!state.serverBusy) {
            state.serverBusy = true;
            double serviceTime = rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, state.clock + serviceTime);
        } else {
            if (params.queueCap == -1 || 
                (int)state.arrivalTimes.size() < params.queueCap) {
                state.arrivalTimes.push(state.clock);
            } else {
                cout << "[ARRIVAL] Customer " << e.customerID 
                     << " rejected (queue full)" << endl;
            }
        }   

        state.numInSystem++;
        cout << "[ARRIVAL] Customer " << e.customerID 
             << " at t=" << fixed << setprecision(4) << e.time << endl;
    }
    
    // ------------------------------------------------------------------------
    // HANDLE DEPARTURE
    // Owner: ANGGOTA 2
    // TODO: Process departure event
    // ------------------------------------------------------------------------
    void handleDeparture(Event e) {
        // TODO ANGGOTA 2:
        // 1. Save previous clock time
        // 2. Update clock to event time
        // 3. Call updateTimeIntegrals(previousTime)
        // 4. Calculate delay for this customer:
        //    - Get arrival time from front of queue (if queue not empty)
        //    - delay = clock - arrivalTime
        // 5. If past warmup: update stats (totalDelay, numServed)
        // 6. numInSystem--
        // 7. Check queue:
        //    a. If queue NOT empty:
        //       - Pop from queue
        //       - Schedule DEPARTURE (time = clock + exponential(mu))
        //    b. If queue empty:
        //       - Set serverBusy = false
        
        cout << "[DEPARTURE] Customer " << e.customerID << " at t=" << e.time << endl;
    }
    
    // ------------------------------------------------------------------------
    // CHECK TERMINATION
    // Owner: ANGGOTA 2
    // TODO: Check if simulation should terminate
    // ------------------------------------------------------------------------
    bool shouldTerminate() {
        // TODO ANGGOTA 2:
        // 1. If BY_SERVED mode: return numServed >= maxServed
        // 2. If BY_TIME mode: return clock >= horizonT
        return false; // PLACEHOLDER
    }
    
    // ------------------------------------------------------------------------
    // CHECK WARMUP
    // Owner: ANGGOTA 2
    // TODO: Check if warmup period has ended
    // ------------------------------------------------------------------------
    bool isWarmupComplete() {
        // TODO ANGGOTA 2:
        // 1. If BY_SERVED: return numServed >= warmup
        // 2. If BY_TIME: return clock >= warmup
        // 3. Mark warmup end time if just completed
        return true; // PLACEHOLDER - assume always warmed up for now
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
        //    a. Get next event from FEL (top + pop)
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
    void scheduleEvent(EventType type, double time) {
        // TODO ANGGOTA 1: Create event and push to FEL
        Event e;
        e.type = type;
        e.time = time;
        e.customerID = nextCustomerID++;
        FEL.push(e);
    }
};

// ============================================================================
// SECTION 5: MULTI-REPLICATION CONTROLLER
// Owner: ANGGOTA 2
// ============================================================================

// ----------------------------------------------------------------------------
// RUN MULTIPLE REPLICATIONS
// TODO ANGGOTA 2: Loop over replications
// ----------------------------------------------------------------------------
vector<RepResult> runReplications(Params params, int numReps) {
    // TODO ANGGOTA 2:
    // 1. Create vector to store results
    // 2. Loop for i = 1 to numReps:
    //    a. Create new Params with seed = baseSeed + i
    //    b. Create DES object
    //    c. Run simulation and get result
    //    d. Store result
    //    e. Print progress
    // 3. Return vector of results
    
    vector<RepResult> results;
    
    for (int i = 1; i <= numReps; i++) {
        cout << "\n=== REPLICATION " << i << " / " << numReps << " ===" << endl;
        
        Params p = params;
        p.seed = params.seed + i;  // Different seed per rep
        
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
    // TODO ANGGOTA 2:
    // 1. Open file for writing
    // 2. Write header: RepID,AvgQ,Utilization,AvgDelay,AvgWait,NumServed,SimTime
    // 3. Write each result as CSV row
    // 4. Close file
    
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
    // TODO ANGGOTA 2:
    // 1. Open file for writing
    // 2. Write header: Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width
    // 3. Write each summary as CSV row
    // 4. Close file
    
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
    // TODO ANGGOTA 2:
    // 1. Create Params with default values
    // 2. Loop through argv:
    //    - Check for --lambda, --mu, --term, --reps, etc
    //    - Parse value after flag
    //    - Update Params
    // 3. Validate parameters (lambda < mu for stability)
    // 4. Return Params
    
    Params p;
    // Default values
    p.lambda = 0.9;
    p.mu = 1.0;
    p.maxServed = 10000;
    p.horizonT = 10000.0;
    p.warmup = 1000;
    p.seed = 12345;
    p.queueCap = -1;  // unlimited
    p.termMode = BY_SERVED;
    p.outdir = "./";
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (arg == "--lambda" && i+1 < argc) {
            p.lambda = stod(argv[++i]);
        }
        else if (arg == "--mu" && i+1 < argc) {
            p.mu = stod(argv[++i]);
        }
        // TODO: Add other parameters
        else if (arg == "--help") {
            printHelp();
            exit(0);
        }
    }
    
    // Validate
    if (p.lambda >= p.mu) {
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
    cout << "  --lambda <value>    Arrival rate (default: 0.9)\n";
    cout << "  --mu <value>        Service rate (default: 1.0)\n";
    // TODO: Add all other options
    cout << "\nExample:\n";
    cout << "  ./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --reps 10\n";
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
    int numReps = 10;  // TODO: Get from CLI
    
    // Print configuration
    cout << "\nConfiguration:" << endl;
    cout << "  Lambda: " << params.lambda << endl;
    cout << "  Mu: " << params.mu << endl;
    cout << "  Rho: " << (params.lambda / params.mu) << endl;
    cout << "  Replications: " << numReps << endl;
    
    // TODO ANGGOTA 2: Run replications
    vector<RepResult> results = runReplications(params, numReps);
    
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