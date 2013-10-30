% This works on 64bit linux
%mex mex_maxflow.cpp ../maxflow.cpp ../graph.cpp 
% Windows, need to give it a compiled boost library
mex -ID:/Code/boost_1_54_0/ -LD:/Code/boost_1_54_0/stage/lib/ -llibboost_system-vc100-mt-1_54 -llibboost_filesystem-vc100-mt-1_54 mex_maxflow.cpp ../maxflow.cpp ../graph.cpp
mex -ID:/Code/boost_1_54_0/ -LD:/Code/boost_1_54_0/stage/lib/ -llibboost_system-vc100-mt-1_54 -llibboost_filesystem-vc100-mt-1_54 mex_maxflow_int.cpp ../maxflow.cpp ../graph.cpp