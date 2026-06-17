/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*
* The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <locale.h>

#include "mode.h"

#if( BUILD_WITH_AMI == 1 && defined( _MSC_VER ) && MSC_AMI == 1 )
#include <malloc.h>
#include <io.h>
#endif


#ifdef _WIN32
#include <crtdbg.h>
#endif
#include "ichimain.h"
#ifdef TARGET_EXE_STANDALONE
#include "inchi_api.h"
#endif

#include "bcf_s.h"
#include "permutation_util.h"

#ifdef RENUMBER_ATOMS_AND_RECALC_V106

/*****************************************************************************/
int rrand(int m)
{
    return
        (int)((double)m * (rand() / (RAND_MAX + 1.0)));
}
/*****************************************************************************/
void shuffle(void* obj, size_t nmemb, size_t size)
{
    void* temp = inchi_malloc(size);
    size_t n = nmemb;
    while (n > 1)
    {
        size_t k = rrand((int)n--);
        if (temp) /* djb-rwth: fixing a NULL pointer dereference */
        {
            memcpy(temp, BYTE(obj) + n * size, size);
            memcpy(BYTE(obj) + n * size, BYTE(obj) + k * size, size);
            memcpy(BYTE(obj) + k * size, temp, size);
        }
    }
#ifdef _WIN32
    _free_dbg(temp, _NORMAL_BLOCK); /* djb-rwth: _free_dbg for _malloc_dbg must be used if Windows SDK is used */
#else
    free(temp); /* djb-rwth: otherwise just free */
#endif
}


/* Use after OrigAtData_Duplicate (permuted <-- saved) */
void OrigAtData_Permute(ORIG_ATOM_DATA* permuted, ORIG_ATOM_DATA* saved, int* numbers)
{
    int i, j, k;
    int nat = saved->num_inp_atoms;
    size_t atsize = sizeof(saved->at[0]);
    for (i = 0; i < nat; i++)
    {
        j = numbers[i];
        memcpy(permuted->at + j, saved->at + i, atsize);
        for (k = 0; k < permuted->at[j].valence; k++)
        {
            permuted->at[j].neighbor[k] = numbers[permuted->at[j].neighbor[k]];
        }
        permuted->at[j].orig_at_number = 1 + numbers[permuted->at[j].orig_at_number - 1];
    }
    if (saved->polymer && permuted->polymer)
    {
        if (saved->polymer->pzz)
        {
            for (k = 0; k < saved->polymer->n_pzz; k++)
            {
                permuted->polymer->pzz[k] = numbers[permuted->polymer->pzz[k]];
            }
        }
        if (saved->polymer->units)
        {
            for (k = 0; k < saved->polymer->n; k++)
            {
                permuted->polymer->units[k]->cap1 = 1 + numbers[permuted->polymer->units[k]->cap1 - 1];
                permuted->polymer->units[k]->cap1 = 1 + numbers[permuted->polymer->units[k]->end_atom1 - 1];
                permuted->polymer->units[k]->cap1 = 1 + numbers[permuted->polymer->units[k]->cap2 - 1];
                permuted->polymer->units[k]->cap1 = 1 + numbers[permuted->polymer->units[k]->end_atom2 - 1];
                if (permuted->polymer->units[k]->alist)
                {
                    for (j = 0; j < permuted->polymer->units[k]->na; j++)
                    {
                        permuted->polymer->units[k]->alist[j] = 1 + numbers[permuted->polymer->units[k]->alist[j] - 1];
                    }
                    for (j = 0; j < permuted->polymer->units[k]->nb; j++)
                    {
                        permuted->polymer->units[k]->blist[2 * j] = 1 + numbers[permuted->polymer->units[k]->blist[2 * j] - 1];
                        permuted->polymer->units[k]->blist[2 * j + 1] = 1 + numbers[permuted->polymer->units[k]->blist[2 * j + 1] - 1];
                    }
                }
            }
        }
    }
    if (saved->v3000 && permuted->v3000)
    {
        if (saved->v3000->atom_index_orig && permuted->v3000->atom_index_orig)
        {
            for (k = 0; k < nat; k++)
            {
                permuted->v3000->atom_index_orig[k] = numbers[permuted->v3000->atom_index_orig[k]];
            }
        }
        if (saved->v3000->atom_index_fin && permuted->v3000->atom_index_fin)
        {
            for (k = 0; k < nat; k++)
            {
                permuted->v3000->atom_index_fin[k] = numbers[permuted->v3000->atom_index_fin[k]];
            }
        }
        if (saved->v3000->n_haptic_bonds && saved->v3000->lists_haptic_bonds && permuted->v3000->n_haptic_bonds && permuted->v3000->lists_haptic_bonds)
        {
            for (j = 0; j < saved->v3000->n_haptic_bonds; j++)
            {
                permuted->v3000->lists_haptic_bonds[j][1] = numbers[permuted->v3000->lists_haptic_bonds[j][1]];
                for (k = 3; k < saved->v3000->lists_haptic_bonds[j][2]; k++)
                {
                    permuted->v3000->lists_haptic_bonds[j][k] = numbers[permuted->v3000->lists_haptic_bonds[j][k]];
                }
            }
        }
        if (saved->v3000->n_steabs && saved->v3000->lists_steabs && permuted->v3000->n_steabs && permuted->v3000->lists_steabs)
        {
            for (j = 0; j < saved->v3000->n_steabs; j++)
            {
                for (k = 2; k < saved->v3000->lists_steabs[j][1] + 2; k++)
                {
                    permuted->v3000->lists_steabs[j][k] = numbers[permuted->v3000->lists_steabs[j][k]];
                }
            }
        }
        if (saved->v3000->n_sterel && saved->v3000->lists_sterel && permuted->v3000->n_sterel && permuted->v3000->lists_sterel)
        {
            for (j = 0; j < saved->v3000->n_sterel; j++)
            {
                for (k = 2; k < saved->v3000->lists_sterel[j][1] + 2; k++)
                {
                    permuted->v3000->lists_sterel[j][k] = numbers[permuted->v3000->lists_sterel[j][k]];
                }
            }
        }
        if (saved->v3000->n_sterac && saved->v3000->lists_sterac && permuted->v3000->n_sterac && permuted->v3000->lists_sterac)
        {
            for (j = 0; j < saved->v3000->n_sterac; j++)
            {
                for (k = 2; k < saved->v3000->lists_sterac[j][1] + 2; k++)
                {
                    permuted->v3000->lists_sterac[j][k] = numbers[permuted->v3000->lists_sterac[j][k]];
                }
            }
        }
    }

    return;
}

#endif