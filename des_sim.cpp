#include <iostream>
#include <queue>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <ctime>
#include <map>
#include <algorithm>
#include <filesystem>
#include <limits>

namespace fs = std::filesystem;

using namespace std;

// ============================================================================
// SECTION 1: ENUMS & CONSTANTS
// ============================================================================
enum EventType { ARRIVAL, DEPARTURE };
enum TerminationMode { BY_SERVED, BY_TIME };

// ============================================================================
// SECTION 2: STRUCTS
// ----------------------------------------------------------------------------
// EVENT STRUCT
struct Event {
    EventType type;
    double time;
    int customerID;  // for logging/trace only

    // reversed for min-heap with priority_queue
    bool operator<(const Event& other) const { return time > other.time; }
};

// ----------------------------------------------------------------------------
// STATE STRUCT
struct State {
    double clock;                  // current simulation time
    int numInSystem;               // queue + service
    bool serverBusy;               // server status
    double nextArrivalTime;        // scheduled next arrival time
    queue<double> arrivalTimes;    // FIFO of arrival times (front = in service)
    queue<int> arrivalIDs;         // parallel queue for IDs

    State() {
        clock = 0.0;
        numInSystem = 0;
        serverBusy = false;
        nextArrivalTime = 0.0;
        while (!arrivalTimes.empty()) arrivalTimes.pop();
        while (!arrivalIDs.empty()) arrivalIDs.pop();
    }
};

// ----------------------------------------------------------------------------
// STATS STRUCT
struct Stats {
    // time-integrals & tallies (post-warmup unless noted)
    double totalDelay;       // sum of time in system W (post-warmup)
    double areaQ;            // ∫ Q(t) dt (queue only, post-warmup)
    double areaB;            // ∫ B(t) dt (busy indicator, post-warmup)
    int    numServed;        // completed customers counted post-warmup
    int    numServedTotal;   // completed customers including warm-up
    int    numArrived;       // accepted arrivals (post-warmup not required)
    double warmupEndTime;    // time when warm-up ended (at a DEPARTURE)

    Stats() {
        totalDelay = 0.0;
        areaQ = 0.0;
        areaB = 0.0;
        numServed = 0;
        numServedTotal = 0;
        numArrived = 0;
        warmupEndTime = 0.0;
    }
};

// ----------------------------------------------------------------------------
// PARAMS STRUCT
struct Params {
    double lambda;              // arrival rate
    double mu;                  // service rate
    int maxServed;              // BY_SERVED mode target (post-warmup!)
    double horizonT;            // BY_TIME mode horizon (post-warmup time)
    int warmup;                 // warm-up in customers completed
    int seed;                   // RNG seed
    int queueCap;               // waiting-room capacity (-1 = unlimited)
    TerminationMode termMode;   // termination mode
    int reps;                   // replications
    string outdir;              // output directory

    Params() {
        lambda = 0.0;
        mu = 1.0;
        maxServed = 10000;
        horizonT = 10000.0;
        warmup = 10;
        seed = static_cast<int>(time(nullptr));
        queueCap = -1;
        termMode = BY_TIME;
        reps = 10;
        outdir = "./";
    }

    bool queueNotExploding() const { return (lambda < mu); }
};

// ----------------------------------------------------------------------------
// REPLICATION RESULT STRUCT
struct RepResult {
    int repID;
    double avgQ;          // Lq
    double utilization;   // ρ
    double avgDelay;      // W
    double avgWait;       // Wq = W - 1/mu
    int numServed;        // post-warmup completed
    double simTime;       // post-warmup time span used for metrics

    RepResult() {
        repID = 0;
        avgQ = 0.0;
        utilization = 0.0;
        avgDelay = 0.0;
        avgWait = 0.0;
        numServed = 0;
        simTime = 0.0;
    }
};

// ============================================================================
// SECTION 3: RANDOM NUMBER GENERATOR CLASS
// ============================================================================
class RNG {
private:
    mt19937_64 generator;

public:
    RNG(int seed) { generator.seed(seed); }

    double exponential(double rate) {
        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator);
        return -log(u) / rate;
    }

    void test(double rate, int samples = 1000) {
        double sum = 0.0;
        for (int i = 0; i < samples; i++) sum += exponential(rate);
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
    priority_queue<Event> FEL;  // Future Event List (min-heap by time)
    int nextCustomerID;

public:
    DES(Params p) : params(p), rng(p.seed) { nextCustomerID = 1; }

    // INIT
    void init() {
        state = State();
        stats = Stats();
        nextCustomerID = 1;
        while (!FEL.empty()) FEL.pop();

        // schedule first arrival
        double firstArrival = rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, firstArrival, nextCustomerID++);
        state.nextArrivalTime = firstArrival;

        cout << "[INIT] t=0, first ARRIVAL at t=" << fixed << setprecision(6)
             << firstArrival << "\n";
    }

    // Helper: are we past warm-up?
    bool isWarmupComplete() const { return stats.numServedTotal >= params.warmup; }

    // UPDATE TIME INTEGRALS (uses state before applying event effects)
    void updateTimeIntegrals(double previousTime) {
        double dt = state.clock - previousTime;
        if (dt <= 0) return;

        if (isWarmupComplete()) {
            int numInQueue = state.numInSystem - (state.serverBusy ? 1 : 0);
            if (numInQueue < 0) numInQueue = 0;
            stats.areaQ += numInQueue * dt;
            stats.areaB += (state.serverBusy ? 1.0 : 0.0) * dt;
        }
    }

    // HANDLE ARRIVAL
    void handleArrival(const Event& e) {
        double prevTime = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(prevTime);

        // Schedule next arrival
        double nextArrival = state.clock + rng.exponential(params.lambda);
        scheduleEvent(ARRIVAL, nextArrival, nextCustomerID++);
        state.nextArrivalTime = nextArrival;

        // Capacity check: only matters if server busy (waiting room)
        bool accepted = true;
        if (state.serverBusy) {
            int waiting = state.numInSystem - 1;
            if (params.queueCap != -1 && waiting >= params.queueCap)
                accepted = false;
        }

        if (!accepted) {
            cout << "[ARRIVAL] Customer " << e.customerID
                 << " rejected at t=" << fixed << setprecision(6) << e.time
                 << " (queue full)\n";
            return; // arrival not admitted, do not change counts
        }

        // Admit: push into queues (front represents customer in service)
        state.arrivalTimes.push(state.clock);
        state.arrivalIDs.push(e.customerID);
        state.numInSystem++;
        if (isWarmupComplete()) stats.numArrived++;

        // If server idle, start service immediately for this (front) customer
        if (!state.serverBusy) {
            state.serverBusy = true;
            double serviceTime = rng.exponential(params.mu);
            scheduleEvent(DEPARTURE, state.clock + serviceTime,
                          state.arrivalIDs.front());
        }

        cout << "[ARRIVAL] Customer " << e.customerID
             << " at t=" << fixed << setprecision(6) << e.time << "\n";
    }

    // HANDLE DEPARTURE
    void handleDeparture(const Event& e) {
        double prevTime = state.clock;
        state.clock = e.time;
        updateTimeIntegrals(prevTime);

        // The customer departing must be the one at the front of the queues
        if (state.arrivalTimes.empty() || state.arrivalIDs.empty()) {
            // Should not happen in a correct M/M/1 flow, guard anyway
            cerr << "[WARN] Departure with empty queue at t=" << e.time << "\n";
        } else {
            double arrived_at = state.arrivalTimes.front();
            int cid = state.arrivalIDs.front();
            (void)cid; // for optional debugging
            state.arrivalTimes.pop();
            state.arrivalIDs.pop();

            double delay = (e.time - arrived_at); // time in system (W)

            // Increase total served (including warm-up)
            stats.numServedTotal++;

            // If this exact departure completes warm-up, mark the time
            if (stats.numServedTotal == params.warmup) {
                stats.warmupEndTime = state.clock; // warm-up ends now
            }

            // Record metrics only post-warmup
            if (isWarmupComplete()) {
                stats.totalDelay += delay;
                stats.numServed++;
            }
        }

        // One customer left the system
        state.numInSystem--;

        // Start next service if anyone waiting
        if (state.numInSystem > 0) {
            // still busy, start service for next in queue
            double serviceTime = rng.exponential(params.mu);
            int nextCID = state.arrivalIDs.front();
            scheduleEvent(DEPARTURE, state.clock + serviceTime, nextCID);
            state.serverBusy = true;
        } else {
            state.serverBusy = false;
        }

        cout << "[DEPARTURE] at t=" << fixed << setprecision(6) << e.time << "\n";
    }

    // TERMINATION CHECK
    bool shouldTerminate() const {
        if (!isWarmupComplete()) return false; // never terminate before warm-up ends
        if (params.termMode == BY_SERVED) {
            return stats.numServed >= params.maxServed;
        } else { // BY_TIME
            double elapsed = state.clock - stats.warmupEndTime;
            return elapsed >= params.horizonT;
        }
    }

    // RUN
    RepResult run() {
        init();

        while (true) {
            if (FEL.empty()) {
                cerr << "[FATAL] FEL became empty. Stopping.\n";
                break;
            }
            if (shouldTerminate()) break;

            Event e = FEL.top(); FEL.pop();
            if (e.type == ARRIVAL)      handleArrival(e);
            else /* DEPARTURE */        handleDeparture(e);
        }

        return computeResults();
    }

    // COMPUTE RESULTS
    RepResult computeResults() {
        RepResult res;
        res.repID = 0;

        double period = max(0.0, state.clock - stats.warmupEndTime);
        if (period > 0.0) {
            res.avgQ = stats.areaQ / period;            // Lq
            res.utilization = stats.areaB / period;     // ρ
        } else {
            res.avgQ = 0.0;
            res.utilization = 0.0;
        }

        if (stats.numServed > 0) {
            res.avgDelay = stats.totalDelay / stats.numServed; // W
        } else {
            res.avgDelay = 0.0;
        }
        res.avgWait = res.avgDelay - (1.0 / params.mu);        // Wq

        res.numServed = stats.numServed;
        res.simTime = period;

        return res;
    }

    // SCHEDULE EVENT
    void scheduleEvent(EventType type, double time, int cid) {
        Event e;
        e.type = type;
        e.time = time;
        e.customerID = cid; // for ARRIVAL we pass a unique ID; for DEPARTURE we pass front's ID
        FEL.push(e);
    }
};

// ============================================================================
// SECTION 5: MULTI-REPLICATION CONTROLLER
// ============================================================================
vector<RepResult> runReplications(Params params, int numReps) {
    vector<RepResult> results;
    results.reserve(numReps);

    for (int i = 1; i <= numReps; i++) {
        cout << "\n=== REPLICATION " << i << " / " << numReps << " ===" << endl;

        Params p = params;
        p.seed = params.seed + i * i;  // distinct seed

        DES sim(p);
        RepResult res = sim.run();
        res.repID = i;
        results.push_back(res);

        cout << "Rep " << i << " complete: "
             << "AvgQ(Lq)=" << res.avgQ
             << " Util(rho)=" << res.utilization
             << " W=" << res.avgDelay
             << " Wq=" << res.avgWait
             << " nServed=" << res.numServed
             << " T=" << res.simTime << endl;
    }
    return results;
}

// ============================================================================
// SECTION 6: STATISTICS & CONFIDENCE INTERVAL
// ============================================================================
struct Summary {
    string metric;
    double mean;
    double stdDev;
    double ci_lower;
    double ci_upper;
    double ci_width;
};

double calculateMean(const vector<double>& data) {
    if (data.empty()) return 0.0;
    double sum = 0.0;
    for (double v : data) sum += v;
    return sum / data.size();
}

double calculateStdDev(const vector<double>& data, double mean) {
    if (data.size() <= 1) return 0.0;
    double sumSq = 0.0;
    for (double v : data) sumSq += (v - mean) * (v - mean);
    double variance = sumSq / (data.size() - 1);
    return sqrt(max(0.0, variance));
}

static double tValue95(int n) {
    // small lookup for 95% two-sided CI
    // df = n-1
    static map<int,double> t = {
        {5,2.776},{6,2.571},{7,2.447},{8,2.365},{9,2.306},
        {10,2.262},{11,2.228},{12,2.201},{13,2.179},{14,2.160},
        {15,2.145},{16,2.131},{17,2.120},{18,2.110},{19,2.101},
        {20,2.093},{21,2.086},{22,2.080},{23,2.074},{24,2.069},
        {25,2.064},{26,2.060},{27,2.056},{28,2.052},{29,2.048},
        {30,2.045}
    };
    if (n >= 30) return 1.96;
    auto it = t.find(n);
    if (it != t.end()) return it->second;
    // for n<5 (degenerate), fall back conservatively
    return 2.776;
}

Summary calculateCI(const vector<double>& data, const string& metricName) {
    Summary s;
    s.metric = metricName;
    int n = (int)data.size();
    s.mean = calculateMean(data);
    s.stdDev = calculateStdDev(data, s.mean);
    if (n <= 1) {
        s.ci_lower = s.ci_upper = s.mean;
        s.ci_width = 0.0;
        return s;
    }
    double tval = tValue95(n);
    double margin = tval * (s.stdDev / sqrt((double)n));
    s.ci_lower = s.mean - margin;
    s.ci_upper = s.mean + margin;
    s.ci_width = 2 * margin;
    return s;
}

vector<Summary> computeSummaries(const vector<RepResult>& results) {
    vector<Summary> summaries;

    vector<double> avgQs, utils, avgDelays, avgWaits;
    avgQs.reserve(results.size());
    utils.reserve(results.size());
    avgDelays.reserve(results.size());
    avgWaits.reserve(results.size());

    for (const auto& r : results) {
        avgQs.push_back(r.avgQ);
        utils.push_back(r.utilization);
        avgDelays.push_back(r.avgDelay);
        avgWaits.push_back(r.avgWait);
    }

    summaries.push_back(calculateCI(avgQs, "AvgQ(Lq)"));
    summaries.push_back(calculateCI(utils, "Utilization(rho)"));
    summaries.push_back(calculateCI(avgDelays, "AvgDelay(W)"));
    summaries.push_back(calculateCI(avgWaits, "AvgWait(Wq)"));

    return summaries;
}

// ============================================================================
// SECTION 7: CSV OUTPUT
// ============================================================================
void ensureOutdirExists(const string& dir) {
    if (dir.empty()) return;
    std::error_code ec;
    fs::create_directories(fs::path(dir), ec); // no-op if it already exists
    if (ec) {
        cerr << "Warning: could not create directory '" << dir
             << "': " << ec.message() << "\n";
    }
}

void writePerRepCSV(const vector<RepResult>& results, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open " << filename << endl;
        return;
    }
    file << "RepID,AvgQ,Utilization,AvgDelay,AvgWait,NumServed,SimTime\n";
    file << fixed << setprecision(8);
    for (const auto& r : results) {
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

void writeSummaryCSV(const vector<Summary>& summaries, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open " << filename << endl;
        return;
    }
    file << "Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width\n";
    file << fixed << setprecision(8);
    for (const auto& s : summaries) {
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

void appendTheoryToSummaryCSV(const Params& p, const string& filename) {
    // ρ is defined even if unstable; other metrics only if λ < μ
    double rho = (p.mu > 0.0) ? (p.lambda / p.mu) : std::numeric_limits<double>::quiet_NaN();

    ofstream file(filename, ios::app);
    if (!file.is_open()) {
        cerr << "Error: Cannot append theory to " << filename << endl;
        return;
    }
    file << fixed << setprecision(8);

    file << "rho(theory)," << rho << ",,,,\n";
    if (p.lambda > 0.0 && p.lambda < p.mu) {
        double W  = 1.0 / (p.mu - p.lambda);
        double Wq = rho / (p.mu - p.lambda);      // = λ / (μ(μ-λ))
        double L  = rho / (1.0 - rho);
        double Lq = (rho * rho) / (1.0 - rho);    // = λ^2 / (μ(μ-λ))
        file << "L(theory),"  << L  << ",,,,\n";
        file << "W(theory),"  << W  << ",,,,\n";
        file << "Wq(theory)," << Wq << ",,,,\n";
        file << "Lq(theory)," << Lq << ",,,,\n";
    } else {
        file << "L(theory),nan,,,,\n";
        file << "W(theory),nan,,,,\n";
        file << "Wq(theory),nan,,,,\n";
        file << "Lq(theory),nan,,,,\n";
    }
    file.close();
}

// ============================================================================
// SECTION 8: CLI PARSER
// ============================================================================
void printHelp();

Params parseArguments(int argc, char* argv[]) {
    Params p;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg == "--help") {
            printHelp();
            exit(0);
        }

        if (i + 1 < argc) {
            if (arg == "--lambda")        p.lambda    = stod(argv[++i]);
            else if (arg == "--mu")       p.mu        = stod(argv[++i]);
            else if (arg == "--horizonT") p.horizonT  = stod(argv[++i]);
            else if (arg == "--maxServed")p.maxServed = stoi(argv[++i]);
            else if (arg == "--warmup")   p.warmup    = stoi(argv[++i]);
            else if (arg == "--reps")     p.reps      = stoi(argv[++i]);
            else if (arg == "--seed")     p.seed      = stoi(argv[++i]);
            else if (arg == "--queueCap") p.queueCap  = stoi(argv[++i]);
            else if (arg == "--term") {
                string t = argv[++i];
                if (t == "served") p.termMode = BY_SERVED;
                else               p.termMode = BY_TIME;
            }
            else if (arg == "--outdir")   p.outdir    = argv[++i];
        }
    }

    if (!p.queueNotExploding()) {
        cerr << "WARNING: Unstable system (lambda >= mu). Ensure this is intended.\n";
    }
    if (p.termMode == BY_SERVED && p.maxServed <= 0) {
        cerr << "WARNING: BY_SERVED requires positive --maxServed.\n";
    }
    if (p.termMode == BY_TIME && p.horizonT <= 0) {
        cerr << "WARNING: BY_TIME requires positive --horizonT.\n";
    }
    return p;
}

void printHelp() {
    cout << "Usage: ./des_sim [OPTIONS]\n\n";
    cout << "Options:\n";
    cout << "  --lambda <v>       Arrival rate λ (default: 0.0)\n";
    cout << "  --mu <v>           Service rate μ (default: 1.0)\n";
    cout << "  --horizonT <v>     Time horizon AFTER warm-up (default: 10000.0)\n";
    cout << "  --maxServed <n>    Customers to complete AFTER warm-up (default: 10000)\n";
    cout << "  --warmup <n>       Customers completed for warm-up deletion (default: 10)\n";
    cout << "  --reps <n>         Number of replications (default: 10)\n";
    cout << "  --seed <n>         RNG seed (default: time-based)\n";
    cout << "  --queueCap <n>     Waiting-room capacity (-1=unlimited) (default: -1)\n";
    cout << "  --term <served|time> Termination mode (default: time)\n";
    cout << "  --outdir <path>    Output directory for CSVs (default: ./)\n\n";
    cout << "Example (BY_SERVED):\n";
    cout << "  ./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --warmup 1000 --reps 10 \\\n";
    cout << "            --seed 12345 --queueCap -1 --term served --outdir ./\n\n";
    cout << "Example (BY_TIME):\n";
    cout << "  ./des_sim --lambda 0.9 --mu 1.0 --horizonT 20000 --warmup 1000 --reps 10 \\\n";
    cout << "            --seed 12345 --queueCap -1 --term time --outdir ./\n";
}

// ============================================================================
// SECTION 9: VALIDATION & ANALYSIS
// ============================================================================
static const Summary* findSummary(const vector<Summary>& S, const string& name) {
    for (const auto& s : S) if (s.metric == name) return &s;
    return nullptr;
}

void validateResults(const Params& p, const vector<Summary>& summaries) {
    double rho = p.lambda / p.mu;
    double L_theory  = rho / (1.0 - rho);
    double W_theory  = 1.0 / (p.mu - p.lambda);
    double Wq_theory = rho / (p.mu - p.lambda);

    const Summary* Lq = findSummary(summaries, "AvgQ(Lq)");
    const Summary* R  = findSummary(summaries, "Utilization(rho)");
    const Summary* W  = findSummary(summaries, "AvgDelay(W)");
    const Summary* Wq = findSummary(summaries, "AvgWait(Wq)");

    cout << "\n=== THEORETICAL VS SIMULATION ===\n" << fixed << setprecision(6);
    cout << "rho (theory)    : " << rho << "\n";
    if (R) cout << "rho (sim mean)  : " << R->mean << "\n";

    // L_sim = Lq + rho
    if (Lq && R) {
        double L_sim = Lq->mean + R->mean;
        cout << "L  (theory)     : " << L_theory << "\n";
        cout << "L  (sim mean)   : " << L_sim
             << "   [via Lq+rho = " << Lq->mean << " + " << R->mean << "]\n";
    } else {
        cout << "L  (theory)     : " << L_theory << "\n";
    }

    cout << "W  (theory)     : " << W_theory << "\n";
    if (W) cout << "W  (sim mean)   : " << W->mean << "\n";

    cout << "Wq (theory)     : " << Wq_theory << "\n";
    if (Wq) cout << "Wq (sim mean)   : " << Wq->mean << "\n";
}

void verifyLittlesLaw(const Params& p, const vector<Summary>& summaries) {
    const Summary* Lq = findSummary(summaries, "AvgQ(Lq)");
    const Summary* R  = findSummary(summaries, "Utilization(rho)");
    const Summary* W  = findSummary(summaries, "AvgDelay(W)");
    if (!(Lq && R && W)) return;

    double L_sim = Lq->mean + R->mean;   // L = Lq + ρ
    double W_sim = W->mean;
    double lambda_used = p.lambda;       // if finite buffer/reneging, use effective rate instead

    double L_from_Little = lambda_used * W_sim;

    cout << "\n=== LITTLE'S LAW VERIFICATION ===\n" << fixed << setprecision(6);
    cout << "L (sim)          = " << L_sim << "\n";
    cout << "λ * W            = " << lambda_used << " * " << W_sim
         << " = " << L_from_Little << "\n";
    cout << "Abs diff         = " << fabs(L_sim - L_from_Little) << "\n";
}

// ============================================================================
// SECTION 10: MAIN
// ============================================================================
int main(int argc, char* argv[]) {
    cout << "==================================================\n";
    cout << "  M/M/1 Queue Discrete Event Simulation\n";
    cout << "==================================================\n";

    Params params = parseArguments(argc, argv);

    cout << "\nConfiguration:\n";
    cout << "  Lambda:       " << params.lambda << "\n";
    cout << "  Mu:           " << params.mu << "\n";
    cout << "  Rho:          " << (params.lambda / params.mu) << "\n";
    cout << "  Termination:  " << (params.termMode == BY_SERVED ? "BY_SERVED" : "BY_TIME") << "\n";
    cout << "  Warm-up (n):  " << params.warmup << "\n";
    cout << "  Replications: " << params.reps << "\n";

    vector<RepResult> results = runReplications(params, params.reps);
    vector<Summary> summaries = computeSummaries(results);

    cout << "\n=== SUMMARY STATISTICS (mean ± CI95) ===\n" << fixed << setprecision(6);
    for (const auto& s : summaries) {
        cout << "  " << s.metric << ": " << s.mean
             << "  [" << s.ci_lower << ", " << s.ci_upper << "]"
             << "  (±" << s.ci_width / 2.0 << ")\n";
    }

    ensureOutdirExists(params.outdir);

    writePerRepCSV(results, params.outdir + "results_per_rep.csv");
    writeSummaryCSV(summaries, params.outdir + "summary.csv");

    appendTheoryToSummaryCSV(params, params.outdir + "summary.csv");

    validateResults(params, summaries);
    verifyLittlesLaw(params, summaries);

    cout << "\n=== SIMULATION COMPLETE ===\n";
    return 0;
}
