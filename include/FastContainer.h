#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

template <typename T>
class RangeStructure
{
private:
  static const int _maxsize = 5;
  int _size;
  std::array<int, _maxsize> _indices;
  std::array<T, _maxsize> _values;
  int _lID = -1;
  int _rID = -1;
public:
  RangeStructure() {
    _size = 0; _indices.fill(-1); _values.fill(std::numeric_limits<T>::quiet_NaN());
  };

  void push_back(const int index, const T& value) {
    int idx = _size++; _indices[idx] = index; _values[idx] = value;
  };

  inline static const int getBatchSize() {return _maxsize;};
  inline const int getSize() const {return _size;}
  inline const int getLNearest() const {return _lID;}
  inline const int getRNearest() const {return _rID;}
  inline void setLNearest(int id) {_lID = id;}
  inline void setRNearest(int id) {_rID = id;}
  inline const std::array<T, _maxsize>& getValues() const {return _values;};
  inline const T getFirst() const {return _values.at(0);};
  inline const T getLast() const {return _values.at(_size-1);};
  inline const std::array<int, _maxsize>& getIndices() const {return _indices;};
  inline const T getFirstID() const {return _indices.at(0);};
  inline const T getLastID() const {return _indices.at(_size-1);};
};

class FastContainer
{
private:
  int _maxSize = RangeStructure<double>::getBatchSize();
  double _lowerBound;
  double _upperBound;
  double _deltaZ = std::numeric_limits<double>::infinity();
  std::vector<RangeStructure<double>> _vec;
public:
    FastContainer(const std::vector<std::pair<int, double>>& input, double lowerBound, double upperBound);
    ~FastContainer() = default;

    // try to make all methods are constant (read safely)
    const int getClosestId(double z) const;
    const std::vector<int> getIdsInRange(double lowerZ, double upperZ) const;
    inline bool isEmpty() const{return !_vec.size();};
};
