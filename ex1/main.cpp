#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>


class ConcludeNext;
class LinearProgramFeasable{
public:
  // dimension of the matrix
  int m, n;
  std::vector< double > b;
  std::vector< double > A;
  LinearProgramFeasable(int m, int n, std::vector< double > b, std::vector< double > A)
  : m(m), n(n), b(b), A(A) {}
public:
  double matrix(int x, int y) const { return A[x + n * y]; }
  
  static LinearProgramFeasable readFromFile (const char *filename);
  LinearProgramFeasable eliminateFirstCoordinate(std::vector<ConcludeNext> &next) const;
  bool trivialFeasable() {
    for (int i = 0; i < m; ++i) {
      if (b[i] < 0)
	return false;
    }
    return true;
  }
};

class ConcludeNext{
public:
  int m, n;
  std::vector<double> A;
  std::vector<double> b;
public:
  double concludeNext(const std::vector<double> &current) {
    double result = -std::numeric_limits<double>::infinity();
    for(int i = 0; i < m; ++i) {
      double cur = b[i];
      for (int j = 0; j < n; ++j) {
	cur += A[j + i * n]*current[j];
      }
      if (cur > result) result = cur;
    }
    if (std::isinf(result))
      // TODO: calculate maximum instead of minimum
      result = -3000;
    return result;
  }
};


std::vector< double > readVectorFromLine(const std::string &line) {
  std::istringstream stream(line);
  double next;
  std::vector<double> result;
  while(!(stream >> next).fail()) {
    result.push_back(next);
  }
  return result;
}

LinearProgramFeasable LinearProgramFeasable::readFromFile(const char* filename) {
  std::fstream input(filename, std::fstream::in);
  int m, n;
  std::string line;
  std::getline(input, line);
  std::istringstream(line) >> m >> n;
  std::getline(input, line);
  std::vector<double> c = readVectorFromLine(line);
  std::getline(input, line);
  std::vector<double> b = readVectorFromLine(line);
  std::vector<double> A;
  for (int i = 0; i < m; ++i) {
    std::getline(input, line);
    auto next = readVectorFromLine(line);
    A.insert(A.end(), next.begin(), next.end());
  }
  return LinearProgramFeasable(m, n, b, A);
}

LinearProgramFeasable LinearProgramFeasable::eliminateFirstCoordinate(std::vector<ConcludeNext> &next) const {
  std::vector<int> linesWithX0GreaterThan0;
  std::vector<int> linesWithX0SmallerThan0;
  std::vector<int> linesWithX0EqualTo0;
  for (int i = 0; i < m; ++i) {
    if (matrix(0, i) > 0)
      linesWithX0GreaterThan0.push_back(i);
    else if (matrix(0, i) < 0)
      linesWithX0SmallerThan0.push_back(i);
    else
      linesWithX0EqualTo0.push_back(i);
  }
  std::vector<double> newA;
  std::vector<double> newB;
  for (auto it1 = linesWithX0GreaterThan0.begin(); it1 != linesWithX0GreaterThan0.end(); ++it1) {
    for (auto it2 = linesWithX0SmallerThan0.begin(); it2 != linesWithX0SmallerThan0.end(); ++it2) {
      for (int z = 1; z < n; ++z) {
	newA.push_back(matrix(z, *it1) / matrix(0, *it1) - matrix(z, *it2) / matrix(0, *it2));
      }
      newB.push_back(b[*it1]/matrix(0, *it1) - b[*it2]/matrix(0, *it2));      
    }
  }
  for (auto it = linesWithX0EqualTo0.begin(); it != linesWithX0EqualTo0.end(); ++it) {
    for (int z = 1; z < n; ++z) {
      newA.push_back(matrix(z, *it));
    }
    newB.push_back(b[*it]);
  }
  ConcludeNext ln;
  ln.n = n - 1;
  ln.m = linesWithX0SmallerThan0.size();
  for (auto it = linesWithX0SmallerThan0.begin(); it != linesWithX0SmallerThan0.end(); ++it) {
    ln.b.push_back(-b[*it]);
    for (int z = 1; z < n; ++z) {
      ln.A.push_back(-matrix(z, *it) / matrix(0, *it));
    }
  }
  next.push_back(ln);
  return LinearProgramFeasable(linesWithX0SmallerThan0.size()*linesWithX0GreaterThan0.size()+linesWithX0EqualTo0.size(), n - 1, newB, newA);
}

int main(int argc, char **argv) {
  auto l1 = LinearProgramFeasable::readFromFile(argv[1]);
  std::vector<ConcludeNext> calculateBack;
  while (l1.n > 0) { l1 = l1.eliminateFirstCoordinate(calculateBack); }
  if (!l1.trivialFeasable()) {
    std::cout << "empty" << std::endl;
  }
  else {
    std::vector<double> result;
    for (auto it = calculateBack.rbegin(); it != calculateBack.rend(); ++it) {
      double next = it->concludeNext(result);
      result.insert(result.begin(), next);
    }
    for (int i = 0; i < result.size(); ++i) {
      std::cout << result[i] << " ";
    }
    std::cout << std::endl;
  }
}
