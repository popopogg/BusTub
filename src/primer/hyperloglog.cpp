//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog.cpp
//
// Identification: src/primer/hyperloglog.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog.h"
#include <algorithm>
#include <cstdint>
#include <mutex>
#include <shared_mutex>
#include "common/util/hash_util.h"

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLog<KeyType>::HyperLogLog(int16_t n_bits) {
  cardinality_ = 0;
  if (n_bits < 0) {
    b_ = 0;
  } else {
    b_ = n_bits;
  }
  m_ = (1 << n_bits);
  buckets_.resize(m_, 0);
}

/**
 * @brief Function that computes binary.
 *
 * @param[in] hash
 * @returns binary of a given hash
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  /** @TODO(student) Implement this function! */
  // return {0};
  //--------------------------------------------------------------------------------------
  return {std::bitset<BITSET_CAPACITY>(hash)};
}

/**
 * @brief Function that computes leading zeros.
 *
 * @param[in] bset - binary values of a given bitset
 * @returns leading zeros of given binary set
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::PositionOfLeftmostOne(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
  /** @TODO(student) Implement this function! */
  // return 0;
  //--------------------------------------------------------------------------------------
  int i = BITSET_CAPACITY - 1 - b_;
  for (; i >= 0; i--) {
    if (bset[i]) {
      return static_cast<uint64_t>(BITSET_CAPACITY - b_ - i);
    }
  }
  return BITSET_CAPACITY - b_ + 1;
}

/**
 * @brief Adds a value into the HyperLogLog.
 *
 * @param[in] val - value that's added into hyperloglog
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::AddElem(KeyType val) -> void {
  /** @TODO(student) Implement this function! */
  //---------------------------------------------------------------------------------------------
  hash_t h = CalculateHash(val);
  std::bitset<BITSET_CAPACITY> bs = ComputeBinary(h);
  uint64_t buckets_no = (bs >> (BITSET_CAPACITY - b_)).to_ullong();
  auto leftmost = static_cast<uint8_t>(PositionOfLeftmostOne(bs));
  std::unique_lock<std::shared_mutex> lk(rwmtx_);
  buckets_[buckets_no] = buckets_[buckets_no] > leftmost ? buckets_[buckets_no] : leftmost;
}

/**
 * @brief Function that computes cardinality.
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeCardinality() -> void {
  /** @TODO(student) Implement this function! */
  std::shared_lock<std::shared_mutex> lk(rwmtx_);
  double sum = 0.0;
  if (m_ == 0) {
    return;
  }

  for (int32_t j = 0; j < m_; ++j) {
    sum += 1.00 / std::pow(2, static_cast<double>(buckets_[j]));
  }
  cardinality_ = static_cast<size_t>(std::floor(CONSTANT * m_ * m_ / sum));
}

template class HyperLogLog<int64_t>;
template class HyperLogLog<std::string>;

}  // namespace bustub
