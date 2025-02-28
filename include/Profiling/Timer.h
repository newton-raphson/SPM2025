//
// Created by maksbh on 4/2/22.
//

#ifndef PROTEUS_TIMER_H
#define PROTEUS_TIMER_H

#include <string>
#include <utility>
#include "DataTypes.h"
#include "talyfem/talyfem.h"

class Timer{
  std::string m_name;
  double start,end;
 public:
  Timer(std::string  name)
  :m_name(std::move(name)){
    start = MPI_Wtime();
  }
  ~Timer(){
    end = MPI_Wtime();
    TALYFEMLIB::PrintInfo("[TIMERS] [ ",m_name," ] : = ",end - start);
  }

};
#endif //PROTEUS_TIMER_H
