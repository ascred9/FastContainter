#include "FastContainer.h"


FastContainer::FastContainer(const std::vector<std::pair<int, double>>& input, double lowerBound, double upperBound):
    _lowerBound(lowerBound),
    _upperBound(upperBound)
{
    // check bounds
    if (upperBound <= lowerBound)
        throw std::invalid_argument("Incorrect upper and lower bounds");

    if (input.size() < _maxSize)
        _deltaZ = upperBound - lowerBound;

    // create ordered set. O(N*lnN)
    // TODO: add epsilon of measurement
    auto comp = [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs)
    {
        return lhs.second < rhs.second;
    };
        
    std::set<std::pair<int, double>, decltype(comp)> tmp_set (input.begin(), input.end(), comp);

    // find the batch step. O(N)
    if (tmp_set.size() <= _maxSize)
        _deltaZ = upperBound - lowerBound;
    else{
        for (auto it = tmp_set.begin(); it != std::prev(tmp_set.end(), _maxSize-1); ++it)
        {
            auto forward_it = std::next(it, _maxSize-1);

            double delta = forward_it->second - it->second;
            if (delta < _deltaZ)
                _deltaZ = delta;
        }
        if (_deltaZ <= 0)
            _deltaZ = upperBound - lowerBound;
    }

    // check bounds against input
    if (tmp_set.begin()->second < lowerBound || std::prev(tmp_set.end())->second > upperBound)
        throw std::invalid_argument("Input is out of range [lower, upper]"); 

    // fill every cell
    _vec = std::vector<RangeStructure<double>>( std::ceil((upperBound - lowerBound) / _deltaZ) );

    // fill new map by input values. O(N)
    for (const auto& elem: tmp_set)
    {
        const int key = (elem.second - lowerBound) / _deltaZ;
        _vec.at(key).push_back(elem.first, elem.second);
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

const int FastContainer::getClosestId(double z) const
{
    int key = (z - _lowerBound) / _deltaZ;

    const auto& p =_vec.at(key);

    if (p.getSize() == 0)
    {
        const int lID = p.getLNearest();
        const int rID = p.getRNearest();
        if (lID > -1 && rID > -1)
            return z - _vec.at(lID).getLast() < _vec.at(rID).getFirst() - z ?  _vec.at(lID).getLastID() : _vec.at(rID).getFirstID();
        else if (lID > -1 && rID == -1)
            return _vec.at(lID).getLastID();
        else if (lID == -1 && rID > -1)
            return _vec.at(rID).getFirstID();
    }

    if (z < p.getFirst())
    {
        const int lID = p.getLNearest();
        return lID > -1 && z - _vec.at(lID).getLast() < p.getFirst() - z ?  _vec.at(lID).getLastID() : p.getFirstID();
    }
    else if (z > p.getLast())
    {
        const int rID = p.getRNearest();
        return rID > -1 && z - p.getLast() > _vec.at(rID).getFirst() - z ? _vec.at(rID).getFirstID() : p.getLastID();
    }

    auto values = p.getValues();
    auto begin = values.begin();
    auto end = values.begin() + p.getSize();
    auto it = std::min_element(begin, end, [&z](const auto& lhs, const auto& rhs){
        return std::isnan(rhs) ? false : std::abs(lhs - z) < std::abs(rhs - z);
    });
    // TODO: lower_bound has O(lnN) for sorted elems.
    auto id = std::distance(values.begin(), it);

    return p.getIndices().at(id);
}

const std::vector<int> FastContainer::getIdsInRange(double lowerZ, double upperZ) const
{
    std::vector<int> res;
    std::vector<double> vals;

    int lowerKey = (lowerZ - _lowerBound) / _deltaZ;
    int upperKey = (upperZ - _lowerBound) / _deltaZ;
    int lastKey = upperKey == _vec.size() ? upperKey : upperKey + 1;

    for (int ikey = lowerKey; ikey < lastKey; ++ikey)
    {
        const auto& p =_vec.at(ikey);
        if (p.getSize() == 0)
            continue;

        res.insert(res.end(), p.getIndices().begin(), p.getIndices().begin() + p.getSize());
        vals.insert(vals.end(), p.getValues().begin(), p.getValues().begin() + p.getSize());
    }

    auto lit = std::lower_bound(vals.begin(), vals.end(), lowerZ);
    auto rit = std::upper_bound(vals.begin(), vals.end(), upperZ);
    auto lDist = std::distance(vals.begin(), lit);
    auto rDist = std::distance(rit, vals.end());
    return {res.begin()+lDist, res.end()-rDist};
}
