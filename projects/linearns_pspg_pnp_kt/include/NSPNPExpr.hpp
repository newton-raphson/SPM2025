//
// Created by skim on 11/21/18.
//
#pragma once

#include <talyfem/talyfem.h>

class NSPNPExpr {
public:
  NSPNPExpr() {
  }

  void set_variable(const double t, const double x, const double y, const double z) {
    Expression::global_symbol_table().set_variable("x", x);
    Expression::global_symbol_table().set_variable("y", y);
    Expression::global_symbol_table().set_variable("z", z);
    Expression::global_symbol_table().set_variable("t", t);
  };

  void set_variable(const double t, const NODE* node) {
    double x = node->location().x();
    double y = node->location().y();
    double z = node->location().z();
    set_variable(t,x,y,z);
  };

  void set_variable(const double t, const FEMElm fe)  {
    double x = fe.position().x();
    double y = fe.position().y();
    double z = fe.position().z();
    set_variable(t,x,y,z);
  }

private:
};
