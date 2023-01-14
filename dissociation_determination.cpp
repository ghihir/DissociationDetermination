#include <bits/stdc++.h>
using namespace std;

struct Atom {
  string element;
  double x, y, z;
};

struct Molecule {
  vector<Atom> atoms;

  double compute_bond_length(int i, int j) {
    double dx = atoms[i].x - atoms[j].x;
    double dy = atoms[i].y - atoms[j].y;
    double dz = atoms[i].z - atoms[j].z;
    return sqrt(dx * dx + dy * dy + dz * dz);
  }

  vector<vector<int>> to_graph() {
    vector<vector<int>> adjacency_list(atoms.size(), vector<int>());
    for (int i = 0; i < atoms.size() - 1; i++) {
      for (int j = i + 1; j < atoms.size(); j++) {
        // auto atom1 = atoms[i], atom2 = atoms[j];
        double length = compute_bond_length(i, j);
        if (is_bonded(length, atoms[i].element, atoms[j].element)) {
          adjacency_list[i].emplace_back(j);
          adjacency_list[j].emplace_back(i);
        }
      }
    }
    return adjacency_list;
  };
};

bool is_bonded(double length, string element1, string element2) {
  map<string, vector<double>> bond_length_list = {
    {"cc", {1.14, 1.62}}, {"ch", {1.02, 1.12}}, {"co", {1.20, 1.50}}, {"cn", {1.40, 1.54}}, {"hh", {0.57, 0.63}}, 
    {"ho", {0.91, 1.01}}, {"hn", {0.95, 1.05}}, {"oo", {1.10, 1.39}}, {"no", {1.14, 1.26}}, {"nn", {1.04, 1.14}}
  };
  string s = element1 + element2;
  sort(s.begin(), s.end());
  return (bond_length_list[s][0] < length && length < bond_length_list[s][0]);
}

struct Molecule {
  vector<string> elements;
  vector<vector<double>> coordinates;

  Molecule (int atom_number) : 
    elements(atom_number), coordinates(atom_number, vector<double>(3)) {}

  double calculate_bond_length(int i, int j) {
    double d_x = coordinates[i][0] - coordinates[j][0];
    double d_y = coordinates[i][1] - coordinates[j][1];
    double d_z = coordinates[i][2] - coordinates[j][2];
    return sqrt(pow(d_x, 2.0) + pow(d_y, 2.0) + pow(d_z, 2.0));
  }

  double get_standard_bond_length(int i, int j) {
    map<string, double> atomic_radius;
    atomic_radius["H"] = 0.32;
    atomic_radius["C"] = 0.75;
    atomic_radius["N"] = 0.71;
    atomic_radius["O"] = 0.63;
    return atomic_radius[elements[i]] + atomic_radius[elements[j]];
  }

  bool is_connected(int i, int j) {
    if (calculate_bond_length(i, j) < get_standard_bond_length(i, j) * 1.3) return true;
    else return false;
  }
};

struct Graph {
  vector<vector<int>> adjacency_list;
  vector<bool> is_visted;
  queue<int> to_visit;
  vector<int> visited_vertexs;

  Graph(int atom_number) : 
    adjacency_list(atom_number), is_visted(atom_number) {}

  void bfs() {
    is_visted[0] = 0;
    to_visit.push(0);
    while (!to_visit.empty()) {
      int current_vertex = to_visit.front();
      to_visit.pop();
      for (auto next_vertex : adjacency_list[current_vertex]) {
        if (is_visted[next_vertex]) continue;
        is_visted[next_vertex] = true;
        to_visit.push(next_vertex);
        visited_vertexs.emplace_back(next_vertex);
      }
    }
  }

  bool is_connected(int atom_number) {
    if (visited_vertexs.size() == atom_number) return true;
    else return false;
  }
};



int main() {
  ifstream fin("input.txt");
  ofstream fout("output.txt");

  int atom_number;
  while (fin >> atom_number) {
    Molecule molecule(atom_number);

    for (int i = 0; i < atom_number; i++) {
      fin >> molecule.elements[i];
      fin >> molecule.coordinates[i][0];
      fin >> molecule.coordinates[i][1];
      fin >> molecule.coordinates[i][2];
    }


    Graph graph(atom_number);
    // cout << graph.adjacency_list.size() << endl;

    for (int i = 0; i < atom_number - 1; i++) {
      for (int j = i + 1; j < atom_number; j++) {
        if (molecule.is_connected(i, j)) {
          graph.adjacency_list[i].emplace_back(j);
          graph.adjacency_list[j].emplace_back(i);
        }
      }
    }

    graph.bfs();

    fout << (graph.is_connected(atom_number) ? 0 : 1) << endl;
  }
}