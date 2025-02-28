#pragma once

#include <stdexcept>

enum NodeDataIndices : int {
  U = 0,
  NODEDATA_MAX = 1
};

class MNodeData {
 public:
  double u_pre[NODEDATA_MAX];

  /**
   * Returns reference to the given value in the object
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  double &value(int index) {
    switch (index) {
      case U: return u_pre[0];
      default: throw std::runtime_error("Invalid MNodeData index");
    }
  }

  /**
   * Const reference version of value().
   * This function is required to be able to get read-only access to values
   * (e.g. when using a `const PPNodeData&` pointer or reference).
   * It is identical to the other value() function except for return type.
   * @param index the index of the desired item
   * @returns const reference to the desired data item
   */
  const double &value(int index) const {
    return const_cast<MNodeData *>(this)->value(index);
  }

  /**
   * Returns the name of the given data value in the object
   * @param index the index of the desired item
   * @return name of the specified data item
   */
  static const char *name(int index) {
    switch (index) {
      case U: return "u";

      default: throw std::runtime_error("Invalid MNodeData index");
    }
  }

  /**
   * @return number of the data items in the object
   */
  static int valueno() {
    return NODEDATA_MAX;
  }
};
