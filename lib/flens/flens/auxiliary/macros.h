/*
 *   Copyright (c) 2009, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FLENS_AUXILIARY_MACROS_H
#define FLENS_AUXILIARY_MACROS_H 1

//-- ADDRESS -------------------------------------------------------------------
#define ADDRESS(x) reinterpret_cast<const void *>(&x)

//-- RAWPOINTER ----------------------------------------------------------------
#define RAWPOINTER(x) reinterpret_cast<const void *>(x)

//-- ASSERT -------------------------------------------------------------------
#include <cassert>

#ifndef ASSERT
#define ASSERT(x) assert(x)
#endif

#ifndef ERROR_MSG
#define ERROR_MSG(x) std::cerr << "ERROR: " << x << std::endl;
#endif

#ifndef NDEBUG

#   ifndef CHECKPOINT_ENTER
#   define CHECKPOINT_ENTER  static bool enter = false; \
                             assert(!enter); enter=true;
#   endif

#   ifndef CHECKPOINT_LEAVE
#   define CHECKPOINT_LEAVE  enter=false;
#   endif

#   ifndef DEBUG_VAR
#   define DEBUG_VAR(x)      x
#   endif

#else

#   ifndef CHECKPOINT_ENTER
#   define CHECKPOINT_ENTER
#   endif

#   ifndef CHECKPOINT_LEAVE
#   define CHECKPOINT_LEAVE
#   endif

#   ifndef DEBUG_VAR
#   define DEBUG_VAR(x)
#   endif

#endif

#endif // FLENS_AUXILIARY_MACROS_H
