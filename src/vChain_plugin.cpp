/*
Copyright 2019 Raffaele Fragapane

http://cultofrig.com
https://github.com/CultOfRig/soft_IK_for_maya

Redistribution and use in source and binary forms, with or without modification, 
  are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer
     in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
     may be used to endorse or promote products
     derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <maya/MFnPlugin.h>

#include "include/vChain_nodes.h"


#if MAYA_API_VERSION < 20180000
    #ifdef _WIN64
    __declspec(dllexport)
    #elif __linux__
    __attribute__((visibility("default")))
    #endif
#endif
MStatus initializePlugin(MObject pluginMob)
{
    MStatus status;
    MFnPlugin fn(pluginMob);

    status = fn.setName("cor_solvers_ik");

    status = fn.registerNode(k_vchain_soft_name, k_vchain_soft_id,
                             &VChain_soft::creator, &VChain_soft::initialize,
                             MPxNode::kDependNode, nullptr);

    return status;
}


#if MAYA_API_VERSION < 20180000
#ifdef _WIN64
    __declspec(dllexport)
#elif __linux__
    __attribute__((visibility("default")))
#endif
#endif
MStatus uninitializePlugin(MObject pluginMob)
{
    MStatus status;


    MFnPlugin fn(pluginMob);

    status = fn.deregisterNode(k_vchain_soft_id);

    return status;
}
