#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

template <typename T>
class FastStructure
{
private:
  static const int _maxsize = 5;
  int _size;
  std::array<std::pair<int, int>, _maxsize> _indices_pos; // first and last index
  std::array<T, _maxsize> _values;
  int _lID = -1;
  int _rID = -1;
public:
  FastStructure() {
    _size = 0; _indices_pos.fill({-1, -1});_values.fill(std::numeric_limits<T>::quiet_NaN());
    };

  void push_back(const int index, const T& value);

  inline static const int getBatchSize() {return _maxsize;};
  inline const int getSize() const {return _size;}
  inline const int getLNearest() const {return _lID;}
  inline const int getRNearest() const {return _rID;}
  inline void setLNearest(int id) {_lID = id;}
  inline void setRNearest(int id) {_rID = id;}
  inline const std::array<T, _maxsize>& getValues() const {return _values;};
  inline const T getFirst() const {return _values.at(0);};
  inline const T getLast() const {return _values.at(_size-1);};
  inline const std::array<std::pair<int, int>, _maxsize>& getIndices() const {return _indices_pos;};
  inline const std::pair<int, int>& getFirstIDpos() const {return _indices_pos.at(0);};
  inline const std::pair<int, int>& getLastIDpos() const {return _indices_pos.at(_size-1);};
};

class FastContainer
{
private:
  int _maxSize = FastStructure<double>::getBatchSize();
  double _lowerBound;
  double _upperBound;
  double _deltaZ = std::numeric_limits<double>::infinity();
  std::vector<FastStructure<double>> _vec;
  std::vector<int> _indices;
public:
    FastContainer() = default;
    FastContainer(double lowerBound, double upperBound);
    ~FastContainer() = default;

    void set(const std::vector<std::pair<int, double>>& input);
    const std::pair<std::vector<int>::const_iterator, std::vector<int>::const_iterator> getClosestId(double z) const;
    const std::pair<std::vector<int>::const_iterator, std::vector<int>::const_iterator> getIdsInRange(double lowerZ, double upperZ) const;
    inline bool isEmpty() const{return !_vec.size();};
};
