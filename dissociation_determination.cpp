#include <bits/stdc++.h>
using namespace std;

struct Atom {
  string element;
  double x, y, z;
};

struct Molecule {
  vector<Atom> atoms;
  vector<vector<int>> graph;
  vector<int> visited;

  double compute_length(Atom atom1, Atom atom2) {
    double dx = atom1.x - atom2.x;
    double dy = atom1.y - atom2.y;
    double dz = atom1.z - atom2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
  }

  bool is_bonded(double length, Atom atom1, Atom atom2) {
    map<string, vector<double>> bond_length_list = {
      {"CC", {1.14, 1.62}}, {"CH", {1.02, 1.12}}, {"CO", {1.20, 1.50}}, {"CN", {1.40, 1.54}}, {"HH", {0.57, 0.63}}, 
      {"HO", {0.91, 1.01}}, {"HN", {0.95, 1.05}}, {"OO", {1.10, 1.39}}, {"NO", {1.14, 1.26}}, {"NN", {1.04, 1.14}}
    };
    transform(atom1.element.begin(), atom1.element.end(), atom1.element.begin(), ::toupper);
    transform(atom2.element.begin(), atom2.element.end(), atom2.element.begin(), ::toupper);
    string s = atom1.element + atom2.element;
    sort(s.begin(), s.end());
    return (bond_length_list[s][0] < length && length < bond_length_list[s][1]);
  }

  void to_graph() {
    graph.resize(atoms.size());
    for (int i = 0; i < atoms.size() - 1; i++) {
      for (int j = i + 1; j < atoms.size(); j++) {
        Atom atom1 = atoms[i], atom2 = atoms[j];
        double length = compute_length(atom1, atom2);
        if (is_bonded(length, atom1, atom2)) {
          graph[i].emplace_back(j);
          graph[j].emplace_back(i);
        }
      }
    }
  };

  void dfs() {
    int s = 0;
    vector<bool> seen(graph.size(), false);
    stack<int> todo;
    visited.emplace_back(s);
    seen[s] = true;
    todo.push(s);
    while (!todo.empty()) {
      int v = todo.top();
      todo.pop();
      for (int x : graph[v]) {
        if (seen[x]) continue;
        visited.emplace_back(x);
        seen[x] = true;
        todo.push(x);
      }
    }
  }

  bool is_dissociated() {
    to_graph();
    dfs();
    return (visited.size() != atoms.size());
  }
};


int main() {
  ifstream fin("input.txt");
  ofstream fout("output.txt");

  int atom_number;
  while (fin >> atom_number) {
    Molecule molecule;

    for (int i = 0; i < atom_number; i++) {
      Atom atom;
      fin >> atom.element >> atom.x >> atom.y >> atom.z;
      molecule.atoms.emplace_back(atom);
    }

    fout << molecule.is_dissociated() << endl;
  }
}