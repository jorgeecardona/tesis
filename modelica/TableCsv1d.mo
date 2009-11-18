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

model Table1D
  Real y;
  Real u;

  parameter Real offset=0; 
  parameter String fileName; // Bug 1116

  function initTable
    input String filename;
    input Real offset;
    output Integer key;
  external "C" 
      annotation(Library="table1d.o", Include="#include \"table1d.h\"");
  end initTable;

  function interpolateTable
    input Integer key;
    input Real u;
    output Real y;
  external "C"
      annotation(Library="table1d.o", Include="#include \"table1d.h\"");
  end interpolateTable;

  function stopTable
    input Integer key;
    output Real y;
  external "C"
      annotation(Library="table1d.o", Include="#include \"table1d.h\"");
  end stopTable;

  Integer key;

equation
  when initial() then
    key = initTable(fileName, offset);
  end when;
  y = interpolateTable(key, u);
end Table1D;

