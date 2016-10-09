#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iterator>
#include <deque>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Graph {
  public:
    using Vertex = size_t;
    using VertexSet = std::unordered_set<Vertex>;
    using AdjencyList = std::unordered_map<Vertex, VertexSet>;

    void AddVertex(Vertex v) {
        adjency_list_[v];
    }

    void AddEdge(Vertex u, Vertex v) {
        adjency_list_[u].insert(v);
        adjency_list_[v].insert(u);
    }

    const VertexSet& AdjecentVertices(Vertex v) const {
        const auto it = adjency_list_.find(v);
        if (it != adjency_list_.end()) {
            return it->second;
        } else {
            return empty_set_;
        }
    }

    VertexSet AllVertices() const {
        VertexSet vs;
        vs.reserve(adjency_list_.size());
        for (const auto& pair : adjency_list_) {
            const auto& vertex = pair.first;
            vs.insert(vertex);
        }
        return vs;
    }

    const AdjencyList& AsAdjencyList() const {
        return adjency_list_;
    }

    AdjencyList& GetAdjencyList() {
        return adjency_list_;
    }

    size_t Size() const {
        return adjency_list_.size();
    }

  private:
    AdjencyList adjency_list_;
    static const VertexSet empty_set_;
};

const Graph::VertexSet Graph::empty_set_;

void GraphEdges(std::ostream& out, const Graph::AdjencyList& adjency_list) {
    for (const auto& pair : adjency_list) {
        const auto& vertex = pair.first;
        const auto& neighbours = pair.second;
        for (const auto adj_vertex : neighbours) {
            out << "\t" << vertex << " -- " << adj_vertex << "\n";
        }
    }
}


struct DebugInfo {
    std::vector<size_t> costs;
};

// Use http://gnuplot.respawned.com/ to plot costs
std::ostream& operator<<(std::ostream& out, const DebugInfo& debug_info) {
    for (size_t i = 0; i < debug_info.costs.size(); ++i) {
        out << i << " " << debug_info.costs[i] << "\n";
    }
    return out;
}

void GoSomePath(Graph::Vertex& v, Graph& graph, std::deque<Graph::Vertex>& path, 
            Graph::VertexSet& comp, std::vector<bool>& used) {
    used[v] = true;
    path.push_back(v);
    comp.erase(v);
    for (auto vertex : graph.GetAdjencyList()[v]) {
        if (!used[vertex])
            GoSomePath(vertex, graph, path, comp, used);
        break;
    }
}

class LongestPath {
 public:
 	explicit LongestPath(Graph& graph)
 		: graph(graph), complement(graph.AllVertices()) {}

    explicit LongestPath(Graph& graph, std::deque<Graph::Vertex>& other_path,
                         Graph::VertexSet& other_comp)
        : graph(graph), path(other_path), complement(other_comp) {}

    void InitRandomPath() {
        auto it = complement.begin();
        for (size_t i = 0; i < rand() % complement.size(); ++i)
            ++it;
        Graph::Vertex v = *it;
        std::vector<bool> used(graph.Size(), false);
        GoSomePath(v, graph, path, complement, used);
    }

    size_t NumAdd() {
        size_t num = 0;
        for (auto elem : graph.GetAdjencyList()[path.front()])
            if (complement.find(elem) != complement.end())
                ++num;
        for (auto elem : graph.GetAdjencyList()[path.back()])
            if (complement.find(elem) != complement.end())
                ++num;
        return num;
    }

    Graph::VertexSet GetNeighbors() {
        if (path.size() == 0) {
            Graph::VertexSet neighbors(complement.begin(), complement.end());
            return neighbors;
        }
        Graph::VertexSet neighbors; 
        neighbors.insert(path.front());
        for (auto elem : graph.GetAdjencyList()[path.front()])
            if (complement.find(elem) != complement.end())
                neighbors.insert(elem);
        if (path.front() != path.back()) {
            neighbors.insert(path.back());
            for (auto elem : graph.GetAdjencyList()[path.back()])
                if (complement.find(elem) != complement.end())
                    neighbors.insert(elem);
        }
    return neighbors;
    }

    const Graph& GetGraph() const {
        return graph;
    }

    void AddFront(const Graph::Vertex& v) {
        path.push_front(v);
        complement.erase(v);
    }
    void AddBack(const Graph::Vertex& v) {
        path.push_back(v);
        complement.erase(v);
    }
    void RemoveFront(const Graph::Vertex& v) {
        path.pop_front();
        complement.insert(v);
    }  
    void RemoveBack(const Graph::Vertex& v) {
        path.pop_back();
        complement.insert(v);
    }

    bool check(const Graph::Vertex& v) {
        for (auto elem : path)
            if (elem == v)
                return true;
        return false;
    }

    void Modify() {
        Graph::VertexSet neighbors = GetNeighbors();
        auto set_it = neighbors.begin();
        for(size_t i = 0; i < rand() % neighbors.size(); ++i)
            ++set_it;
        Graph::Vertex v = *set_it;
        if (path.size() == 0) {
            path.push_back(v);
        } 
        if (v == path.front()) {
            RemoveFront(v);
            return;
        }
        if (v == path.back()) {
            RemoveBack(v);
            return;
        }
        auto FrontAdjList = graph.GetAdjencyList()[path.front()];
        auto BackAdjList = graph.GetAdjencyList()[path.back()];
        if (BackAdjList.find(v) != BackAdjList.end()) {
            AddBack(v);
            return;
        }
        if (FrontAdjList.find(v) != FrontAdjList.end()) {
            AddFront(v);
            return;
        }
    }

    std::deque<Graph::Vertex>& GetPath() {
        return path;
    }

    Graph::VertexSet& GetComp() {
        return complement;
    }

    size_t Cost() const {
        return path.size();
    }

    void PrintPath() {
        for (auto elem : path)
            std::cout << elem << " ";
        std::cout << "\n";
    }
 private:
 	Graph& graph;
 	Graph::VertexSet complement;
    std::deque<Graph::Vertex> path;
};

class LongestPathSolver {
 public:
     virtual LongestPath Solve(Graph& graph,
                               DebugInfo& debug_info) const = 0;
     virtual ~LongestPathSolver() = default;
};

Graph RandomGraph(size_t size, double edge_probability) {
    Graph graph;
    for (Graph::Vertex v = 1; v <= size; ++v) {
        graph.AddVertex(v);
    }
    for (Graph::Vertex v = 1; v <= size; ++v) {
        for (Graph::Vertex u = v + 1; u <= size; ++u) {
            if (double(rand()) / RAND_MAX <= edge_probability) {
                graph.AddEdge(v, u);
            }
        }
    }
    return graph;
}

class GradientDescent final: public LongestPathSolver {
    LongestPath Solve(Graph& graph, DebugInfo&debug_info) const {
        LongestPath lp = LongestPath(graph);
        lp.InitRandomPath();
        debug_info.costs.push_back(lp.Cost());
        while (true) {
            if (lp.NumAdd() == 0)
                break;
            LongestPath next = LongestPath(graph, lp.GetPath(), lp.GetComp());
            next.Modify();
            if (next.Cost() > lp.Cost()) {
                lp.GetPath() = next.GetPath();
                lp.GetComp() = next.GetComp();
            }
            debug_info.costs.push_back(lp.Cost());
        }
        return lp; 
    }
};

class Metropolis final: public LongestPathSolver {
 public:
    Metropolis (double k, double T, bool to_burn)
        : k(k), T(T), to_burn(to_burn) {}
    LongestPath Solve(Graph& graph, DebugInfo&debug_info) const {
        LongestPath lp = LongestPath(graph);
        lp.InitRandomPath();
        debug_info.costs.push_back(lp.Cost());
        for (size_t i = 0; i < 500; ++i) {
            LongestPath next = LongestPath(graph, lp.GetPath(), lp.GetComp());
            next.Modify();
            if (next.Cost() >= lp.Cost()) {
                lp.GetPath() = next.GetPath();
                lp.GetComp() = next.GetComp();               
            } else {
                double rand_num = (double)(rand())/RAND_MAX;
                double E = (double)next.Cost() - lp.Cost();
                double GB = exp(E / (k * T));
                if (rand_num <= GB) {
                    lp.GetPath() = next.GetPath();
                    lp.GetComp() = next.GetComp();                    
                }
            }
            debug_info.costs.push_back(lp.Cost());
            if (to_burn)
                T = T / 1.01;
        }
    return lp;
    }
 private:
    double k;
    mutable double T;
    bool to_burn;
};

int InitRandSeed(int argc, const char* argv[]) {
    int rand_seed;
    if (argc >= 2) {
        rand_seed = atoi(argv[1]);
    } else {
        rand_seed = time(nullptr);
    }
    srand(rand_seed);
    return rand_seed;
}

// Use http://www.webgraphviz.com to take a look at the graph
void GraphViz(std::ostream& out, const Graph& graph) {
    out << "strict graph {\n";
    for (const auto& pair : graph.AsAdjencyList()) {
        const auto& vertex = pair.first;
        out << "\t" << vertex << "\n";
    }
    GraphEdges(out, graph.AsAdjencyList());
    out << "}\n";
}

void GraphViz(std::ostream& out, LongestPath& longest_path) {
    out << "strict graph {\n";
    for (const auto& pair : longest_path.GetGraph().AsAdjencyList()) {
        const auto& vertex = pair.first;
        if (longest_path.check(vertex)) {
            out << "\t" <<  vertex << " [shape=doublecircle]\n";
        } else {
            out << "\t" << vertex << "\n";
        }
    }
    GraphEdges(out, longest_path.GetGraph().AsAdjencyList());
    out << "}\n";
}

void TrySolver(const LongestPathSolver& solver, Graph& graph) {
    GraphViz(std::cout, graph);
    auto best_cost = 0;
    size_t results = 0;
    for (int attempt = 1; attempt < 3; ++attempt) {
        DebugInfo debug_info;
        auto longest_path = solver.Solve(graph, debug_info);
        auto cost = longest_path.Cost();
        if (cost > best_cost) {
            best_cost = cost;
            GraphViz(std::cout, longest_path);
            std::cout << "Trace info:\n" << debug_info << "\n";
            ++results;
        }
    }
    std::cout << "Results: " << results << std::endl;
}

int main(int argc, const char* argv[]) {
    std::cout << "Using rand seed: " << InitRandSeed(argc, argv) << "\n";

    auto graph = RandomGraph(100, 0.03);
    GradientDescent gradient_descent;
    Metropolis metropolis(1, 100, false);
    Metropolis burn(1, 100, true);
    std::cout << "GD:----------------------------------------------\n";
    TrySolver(gradient_descent, graph);
    std::cout << "Metropolis:--------------------------------------\n";
    TrySolver(metropolis, graph);
    std::cout << "Metropolis with burning:-------------------------\n";
    TrySolver(burn, graph);
    return 0;
}

