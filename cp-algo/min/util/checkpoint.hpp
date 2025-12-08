#ifndef CP_ALGO_UTIL_CHECKPOINT_HPP
#define CP_ALGO_UTIL_CHECKPOINT_HPP
#include "../util/big_alloc.hpp"
#include <iostream>
#include <chrono>
#include <string>
#include <map>
namespace cp_algo{
#ifdef CP_ALGO_CHECKPOINT
big_map<big_string,double>checkpoints;double last;
#endif
template<bool final=false>void checkpoint([[maybe_unused]]auto const&_msg){
#ifdef CP_ALGO_CHECKPOINT
big_string msg=_msg;double now=(double)clock()/CLOCKS_PER_SEC;double delta=now-last;last=now;if(msg.size()&&!final){checkpoints[msg]+=delta;}if(final){for(auto const&[key,value]:checkpoints){std::cerr<<key<<": "<<value*1000<<" ms\n";}std::cerr<<"Total: "<<now*1000<<" ms\n";}
#endif
}template<bool final=false>void checkpoint(){checkpoint<final>("");}}
#endif
