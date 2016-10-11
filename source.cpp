#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iterator>
#include <list>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template<typename T>
std::unordered_set<T> SetIntersection (const std::unordered_set<T>& a,
                                       const std::unordered_set<T>& b) {
    std::unordered_set<T> c;
    for (const auto& elem : a)
        if (b.count(elem) != 0)
            c.insert(elem);
    return c;
}

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

std::ostream& operator<<(std::ostream& out, const DebugInfo& debug_info) {
    for (size_t i = 0; i < debug_info.costs.size(); ++i) {
        out << i << " " << debug_info.costs[i] << "\n";
    }
    return out;
}

void GoSomePath(Graph::Vertex& v, const Graph& graph, std::list<Graph::Vertex>& path, 
            Graph::VertexSet& comp, std::vector<bool>& used) {
    used[v] = true;
    path.push_back(v);
    comp.erase(v);
    for (auto vertex : graph.AsAdjencyList().at(v)) {
        if (!used[vertex])
            GoSomePath(vertex, graph, path, comp, used);
        break;
    }
}

class LongestPath {
 public:
 	explicit LongestPath(const Graph& graph)
 		: graph(graph), complement(graph.AllVertices()) {}

    explicit LongestPath(const Graph& graph, std::list<Graph::Vertex>& other_path,
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
        auto path_it_1 = path.begin();
        auto path_it_2 = path.begin();
        
        if (path_it_1 == path.end())
            return complement.size();

        num += SetIntersection(graph.AsAdjencyList().at(*path_it_1), complement).size();

        ++path_it_2;
        if (path_it_2 == path.end())
            return num;
        
        for (; path_it_2 != path.end(); ++path_it_2, ++path_it_1) {
            auto inter = SetIntersection(graph.AsAdjencyList().at(*path_it_1),
                                         graph.AsAdjencyList().at(*path_it_2));
            inter = SetIntersection(inter, complement);
        } 

        num += SetIntersection(graph.AsAdjencyList().at(*path_it_1), complement).size();
        return num;
    }

    const Graph& GetGraph() const {
        return graph;
    }

    bool check(const Graph::Vertex& v) const {
        for (auto elem : path)
            if (elem == v)
                return true;
        return false;
    }

    void Modify() {
        size_t path_size = path.size();
        if (path_size == 0) {
            auto comp_it = complement.begin();
            for (size_t i = 0; i < rand() % complement.size(); ++i)
                ++comp_it;
            complement.erase(*comp_it);
            path.push_back(*comp_it);
            return;
        }
        auto path_it = path.begin();
        for (size_t i = 0; i < rand() % path_size; ++i)
            ++path_it;
        const Graph::Vertex v = *path_it;
        Graph::VertexSet v_adj;
        v_adj = SetIntersection(graph.AsAdjencyList().at(v), complement);
        size_t comp_size = complement.size();
// Рассмотрим случай, когда в пути всего одна вершина.
        if (path_size == 1) {
            // Число соседей - complement_size + v_adj_size + 1
            // т.к. для каждой вершины из complement мы можем заменить v на нее,
            // либо соединить v с каждой вершиной из v_adj, 
            // и один вариант, когда удаляем v
            size_t v_adj_size = v_adj.size();
            size_t NumNeighbors = comp_size + v_adj_size + 1;
            size_t variant = 0;
            if (NumNeighbors != 1)
                variant = rand() % (NumNeighbors - 1); // учтем 0
            if (0 <= variant && variant < comp_size) {
                auto comp_it = complement.begin();
                for (size_t i = 0; i < variant; ++i)
                    ++comp_it;
                complement.erase(*comp_it);
                complement.insert(v);
                path.erase(path_it);
                path.push_back(*comp_it);
                return;
            }
            if (comp_size <= variant && variant <  v_adj_size + comp_size) {
                variant -= comp_size;
                auto v_adj_it = v_adj.begin();
                for(size_t i = 0; i < variant; ++i)
                    ++v_adj_it;
                complement.erase(*v_adj_it);
                path.push_back(*v_adj_it);
                return;
            }
            if (variant == comp_size + v_adj_size) {
                complement.insert(v);
                path.erase(path_it);
                return;
            }
        }
// Теперь случай, когда мы попали в крайнюю вершину
        if (v == path.front() || v == path.back()) {
            bool is_front;
            Graph::Vertex next;
            auto next_it = path_it;
            if (v == path.front()) {
                next = *(++next_it);
                is_front = true;
            }
            if (v == path.back()) {
                next = *(--next_it);
                is_front = false;
            }
            Graph::VertexSet next_adj;
            next_adj = SetIntersection(graph.AsAdjencyList().at(next), complement);

            size_t v_adj_size = v_adj.size();
            size_t next_adj_size = next_adj.size();
    // Количество соседей: next_adj_size + v_adj_size + 1
    // next_adj_size вариантов, чтобы заменить вершину v
    // v_adj_size вариантов для добавления вершины к v
    // 1 вариант, чтобы удалить v
            size_t NumNeighbors = next_adj_size + v_adj_size + 1;
            size_t variant = 0;
            if (NumNeighbors != 1)
                variant = rand() % (NumNeighbors - 1); // учтем 0
//            std::cout << "BEFORE:\n\t";
//            PrintPath();
            if (0 <= variant && variant < next_adj_size) {
                auto next_adj_it = next_adj.begin();
                for (size_t i = 0; i < variant; ++i)
                    ++next_adj_it;
                complement.erase(*next_adj_it);
                complement.insert(v);
                path.erase(path_it);
                if (is_front)
                    path.push_front(*next_adj_it);
                else
                    path.push_back(*next_adj_it);
//                std::cout << v << " заменили на " << *next_adj_it << "\n";
//                std::cout << "AFTER:\n\t";
                PrintPath();
                return;
            }
            if (next_adj_size <= variant && variant <  v_adj_size + next_adj_size) {
                variant -= next_adj_size;
                auto v_adj_it = v_adj.begin();
                for(size_t i = 0; i < variant; ++i)
                    ++v_adj_it;
                complement.erase(*v_adj_it);
                if (is_front)
                    path.push_front(*v_adj_it);
                else
                    path.push_back(*v_adj_it);
//                std::cout << "к " << v << " добавили " << *v_adj_it << "\n";
//                std::cout << "AFTER:\n\t";
                PrintPath();
                return;
            }
            if (variant == next_adj_size + v_adj_size) {
                complement.insert(v);
                path.erase(path_it);
//                std::cout << v << " удалили!\n";
//                std::cout << "AFTER:\n\t";
                PrintPath();
                return;
            }
        }
// Случай, когда мы попали в середину пути
        auto next_it = path_it;
        auto prev_it = path_it;
        auto next = *(++next_it);
        auto prev = *(--prev_it);
        
        Graph::VertexSet next_adj;
        Graph::VertexSet prev_adj;
        next_adj = SetIntersection(graph.AsAdjencyList().at(next), complement);
        prev_adj = SetIntersection(graph.AsAdjencyList().at(prev), complement);

        Graph::VertexSet substitude; // множество вершин, на которые можно заменить v
        substitude = SetIntersection(next_adj, prev_adj);

        Graph::VertexSet prev_add; // множество вершин, которых можно вставить между v и предыдущей вершиной
        prev_add = SetIntersection(v_adj, prev_adj);
       
        Graph::VertexSet next_add; // множество вершин, которых можно вставить между v и следующей вершиной
        next_add = SetIntersection(next_adj, v_adj);

        size_t subs_size = substitude.size();
        size_t next_add_size = next_add.size();
        size_t prev_add_size = prev_add.size();
        bool del_size = 0;
        if (next_adj.count(prev) != 0)
            del_size = 1;

        // получили subs_size + next_add_size + prev_add_size + del_size соседей
        size_t NumNeighbors = subs_size + next_add_size + prev_add_size + del_size;
        size_t variant = 0;
        if (NumNeighbors != 1)
            variant = rand() % (NumNeighbors - 1);
        if (0 <= variant && variant < subs_size) {
            auto sub_it = substitude.begin();
            for (size_t i = 0; i < variant; ++i)
                ++sub_it;
            complement.erase(*sub_it);
            complement.insert(v);
            path.insert(path_it, *sub_it);
            path.erase(path_it);
            return;
        }
        if (subs_size <= variant && variant < subs_size + next_add_size) {
            variant -= subs_size;
            auto next_add_it = next_add.begin();
            for (size_t i = 0; i < variant; ++i)
                ++next_add_it;
            complement.erase(*next_add_it);
            path.insert(next_it, *next_add_it);
            return;
        }
        if (subs_size + next_add_size <= variant && 
                variant < subs_size + next_add_size + prev_add_size) {
            variant -= (subs_size + next_add_size);
            auto prev_add_it = prev_add.begin();
            for (size_t i = 0; i < variant; ++i)
                ++prev_add_it;
            complement.erase(*prev_add_it);
            path.insert(path_it, *prev_add_it);
            return;
        }
        if (subs_size + next_add_size + prev_add_size <= variant &&
                variant <= subs_size + next_add_size + prev_add_size + del_size) {
            path.erase(path_it);
            complement.insert(v);
            return;
        }
    }

    std::list<Graph::Vertex>& GetPath() {
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
 	const Graph& graph;
 	Graph::VertexSet complement;
    std::list<Graph::Vertex> path;
};

class LongestPathSolver {
 public:
     virtual LongestPath Solve(const Graph& graph,
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
    LongestPath Solve(const Graph& graph, DebugInfo&debug_info) const {
        LongestPath lp = LongestPath(graph);
        lp.InitRandomPath();
        debug_info.costs.push_back(lp.Cost());
        while (true) {
            auto num = lp.NumAdd();
            if (num  == 0)
                break;
            LongestPath next = LongestPath(graph, lp.GetPath(), lp.GetComp());
            next.Modify();
            if (next.Cost() > lp.Cost()) {
                lp.GetPath() = next.GetPath();
                lp.GetComp() = next.GetComp();
            }
            debug_info.costs.push_back(lp.Cost());
        }
//        lp.PrintPath();
        return lp; 
    }
};

class Metropolis final: public LongestPathSolver {
 public:
    Metropolis (double k, double T, bool annealing)
        : k(k), T(T), annealing(annealing) {}
    LongestPath Solve(const Graph& graph, DebugInfo&debug_info) const {
        LongestPath lp = LongestPath(graph);
        lp.InitRandomPath();
        debug_info.costs.push_back(lp.Cost());
        for (size_t i = 0; i < 100; ++i) {
            LongestPath next = LongestPath(graph, lp.GetPath(), lp.GetComp());
            next.Modify();
            if (next.Cost() > lp.Cost()) {
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
//            lp.PrintPath();
            debug_info.costs.push_back(lp.Cost());
            if (annealing)
                T = T / 1.023;
        }
//    lp.PrintPath();
    return lp;
    }
 private:
    double k;
    mutable double T;
    bool annealing;
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

void GraphViz(std::ostream& out, const Graph& graph) {
    out << "strict graph {\n";
    for (const auto& pair : graph.AsAdjencyList()) {
        const auto& vertex = pair.first;
        out << "\t" << vertex << "\n";
    }
    GraphEdges(out, graph.AsAdjencyList());
    out << "}\n";
}

void GraphViz(std::ostream& out, const LongestPath& longest_path) {
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

void TrySolver(const LongestPathSolver& solver, const Graph& graph) {
//    GraphViz(std::cout, graph);
    auto best_cost = 0;
    size_t results = 0;
    for (int attempt = 1; attempt < 2; ++attempt) {
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

    auto graph = RandomGraph(10, 0.2);
    GradientDescent gradient_descent;
    Metropolis metropolis(1, 100, false);
    Metropolis ann_simulation(1, 100, true);
    std::cout << "GD:----------------------------------------------\n";
    TrySolver(gradient_descent, graph);
    std::cout << "Metropolis:--------------------------------------\n";
    TrySolver(metropolis, graph);
    std::cout << "Metropolis with annealing:-----------------------\n";
    TrySolver(ann_simulation, graph);
    return 0;
}

