/* 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or    
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of    
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    
GNU General Public License for more details.

You should have received a copy of the GNU General Public License    
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TABLE1D_H
#define TABLE1D_H
// Function that initiliaze the table lookup
int initTable(const char* filename, const double offset);

// Function that inerpolate the table lookup
double interpolateTable(const int key, const double u);

// Function that delete the shared memory
void stopTable(const int key);

#endif
