#include "FastContainer.h"


FastContainer::FastContainer(double lowerBound, double upperBound):
    _lowerBound(lowerBound),
    _upperBound(upperBound)
{
    // check bounds
    if (upperBound <= lowerBound)
        throw std::invalid_argument("Incorrect upper and lower bounds");
}

void FastContainer::set(const std::vector<std::pair<int, double>>& input)
{
    _vec.clear();
    _indices.clear();

    if (input.size() < _maxSize)
        _deltaZ = _upperBound - _lowerBound;

    // create ordered set. O(N*lnN)
    auto comp = [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs)
    {
        return lhs.second < rhs.second;
    };
    std::multiset<std::pair<int, double>, decltype(comp)> tmp_multiset(input.begin(), input.end(), comp);

    // create value set. O(N) because tmp_multiset is already sorted
    std::set<std::pair<int, double>, decltype(comp)> tmp_val_set(tmp_multiset.begin(), tmp_multiset.end(), comp);

    // find the batch step. O(N)
    if (tmp_val_set.size() <= _maxSize)
        _deltaZ = _upperBound - _lowerBound;
    else{
        for (auto it = tmp_val_set.begin(); it != std::prev(tmp_val_set.end(), _maxSize-1); ++it)
        {
            auto forward_it = std::next(it, _maxSize-1);

            double delta = forward_it->second - it->second;
            if (delta < _deltaZ)
                _deltaZ = delta;
        }
        if (_deltaZ <= 0)
            _deltaZ = _upperBound - _lowerBound;
    }

    // check bounds against input
    if (tmp_val_set.begin()->second < _lowerBound || std::prev(tmp_val_set.end())->second > _upperBound)
        throw std::invalid_argument("Input is out of range [lower, upper]"); 

    // fill every cell
    _vec = std::vector<FastStructure<double>>( std::ceil((_upperBound - _lowerBound) / _deltaZ) );
    _indices.reserve(input.size());

    // fill new map by input values. O(N)
    int last_index = 0;
    for (const auto& elem: tmp_multiset)
    {
        const int key = (elem.second - _lowerBound) / _deltaZ;
        _vec.at(key).push_back(last_index, elem.second);
        _indices.emplace_back(elem.first);
        ++last_index;
    }

    // fill empty cells by indices to the nearest
    int tmpL = -1;
    for (int idx = 0; idx < _vec.size(); ++idx)
    {
        _vec.at(idx).setLNearest(tmpL);
        if (_vec.at(idx).getSize() != 0)
            tmpL = idx;
    }

    int tmpR = -1;
    for (int idx = _vec.size() - 1; idx >= 0; --idx)
    {
        _vec.at(idx).setRNearest(tmpR);
        if (_vec.at(idx).getSize() != 0)
            tmpR = idx;
    }
}

const std::pair<std::vector<int>::const_iterator, std::vector<int>::const_iterator> FastContainer::getClosestId(double z) const
{
    int key = (z - _lowerBound) / _deltaZ;

    const auto& p =_vec.at(key);

    if (p.getSize() == 0)
    {
        const int lID = p.getLNearest();
        const int rID = p.getRNearest();
        if (lID > -1 && rID > -1)
        {
            auto pos = z - _vec.at(lID).getLast() < _vec.at(rID).getFirst() - z ?  _vec.at(lID).getLastIDpos() : _vec.at(rID).getFirstIDpos();
            return {_indices.begin() + pos.first, std::next(_indices.begin() + pos.second)};
        }
        else if (lID > -1 && rID == -1)
        {
            auto pos = _vec.at(lID).getLastIDpos();
            return {_indices.begin() + pos.first, std::next(_indices.begin() + pos.second)};
        }
        else if (lID == -1 && rID > -1)
        {
            auto pos = _vec.at(rID).getFirstIDpos();
            return {_indices.begin() + pos.first, std::next(_indices.begin() + pos.second)};
        }
    }

    if (z < p.getFirst())
    {
        const int lID = p.getLNearest();
        auto pos = lID > -1 && z - _vec.at(lID).getLast() < p.getFirst() - z ?  _vec.at(lID).getLastIDpos() : p.getFirstIDpos();
        return {_indices.begin() + pos.first, std::next(_indices.begin() + pos.second)};
    }
    else if (z > p.getLast())
    {
        const int rID = p.getRNearest();
        auto pos = rID > -1 && z - p.getLast() > _vec.at(rID).getFirst() - z ? _vec.at(rID).getFirstIDpos() : p.getLastIDpos();
        return {_indices.begin() + pos.first, std::next(_indices.begin() + pos.second)};
    }

    auto values = p.getValues();
    auto begin = values.begin();
    auto end = values.begin() + p.getSize();
    auto it = std::min_element(begin, end, [&z](const auto& lhs, const auto& rhs){
        return std::isnan(rhs) ? false : std::abs(lhs - z) < std::abs(rhs - z);
    });
    // TODO: lower_bound has O(lnN) for sorted elems.
    auto id = std::distance(values.begin(), it);

    auto pos = p.getIndices().at(id);
    return {_indices.begin() + pos.first, std::next(_indices.begin() + pos.second)};
}

const std::pair<std::vector<int>::const_iterator, std::vector<int>::const_iterator> FastContainer::getIdsInRange(double lowerZ, double upperZ) const
{
    int lowerKey = (lowerZ - _lowerBound) / _deltaZ;
    int upperKey = (upperZ - _lowerBound) / _deltaZ;

    lowerKey = lowerKey < _vec.size() ? lowerKey : _vec.size() - 1;
    lowerKey = lowerKey > -1 ? lowerKey : 0;
    
    upperKey = upperKey < _vec.size() ? upperKey : _vec.size() - 1;
    upperKey = upperKey > -1 ? upperKey : 0;

    int first = 0;
    const auto& pLower = _vec.at(lowerKey);
    if (pLower.getSize() == 0)
    {
        if (pLower.getRNearest() == -1)
            first = -1;
        else
        {
            const auto& pR = _vec.at(pLower.getRNearest());
            first = pR.getSize() > 0 ? pR.getFirstIDpos().first : -1;
        }
    }
    else
    {
        auto lit = pLower.getValues().begin();
        lit = std::lower_bound(pLower.getValues().begin(), pLower.getValues().begin()+pLower.getSize(), lowerZ);
        auto lDist = std::distance(pLower.getValues().begin(), lit);
        first = (pLower.getIndices().begin() + lDist)->first;
    }

    int last = 0;
    const auto& pUpper = _vec.at(upperKey);
    if (pLower.getSize() == 0)
    {
        if (pUpper.getLNearest() == -1)
            last = -1;
        else
        {
            const auto& pL = _vec.at(pUpper.getLNearest());
            last = pL.getSize() > 0 ? pL.getLastIDpos().second : -1;
        }
    }
    else
    {
        auto rit = pUpper.getValues().end();
        rit = std::upper_bound(pUpper.getValues().begin(), pUpper.getValues().begin()+pUpper.getSize(), upperZ);
        auto rDist = std::distance(rit, pUpper.getValues().end());
        last = (pUpper.getIndices().end() - rDist)->second;
    }

    if (last >= first && first != -1 && last != -1)
        return {_indices.begin() + first, std::next(_indices.begin() + last)};
    else
        return {_indices.end(), _indices.end()};
}

template <typename T>
void FastStructure<T>::push_back(const int index, const T& value)
{
    int idx = -1;
    for (int id = 0; id < _size; ++id)
    {
        if (value == _values[id])
            idx = id;
    }
        
    if (idx == -1)
    {
        idx = _size++;
        _values[idx] = value;
        _indices_pos[idx] = {index, index};
    }
    else
    {
        int key = _indices_pos[idx].second + 1;
        _indices_pos[idx].second = index;
    
        if (key != _indices_pos[idx].second)
            throw std::invalid_argument("Error in filling");
    }
}
