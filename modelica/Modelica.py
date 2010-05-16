# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or    
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of    
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License    
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pexpect
from math import sqrt
from random import random
from tempfile import NamedTemporaryFile
from Cheetah.Template import Template
import scipy.optimize

class Modelica:
    def __init__(self):
        print "Initializing Modelica Shell"
        for i in xrange(10):
            try:
                self.shell = pexpect.spawn("OMShell-terminal")
                self.shell.expect(">>>")
            except:
                pass
        self.parameters = {}
        self.opt_parameters = {}
        print "Modelica Shell Initialized"

    def __send__(self, cmd, timeout = -1):
        self.shell.sendline(cmd)
        self.shell.expect(">>>", timeout = timeout)
        return self.shell.before

    def loadModel(self, model):
        print "Loading %s"%(model)
        print self.__send__("loadModel(%s)"%(model))

    def loadFile(self, name):
        print "Loading %s"%(name)
        print self.__send__("loadFile(\"%s\")"%(name))

        
    def listVariables(self):
        print self.__send__("listVariables()")

    def base_simulate(self, model, stopTime = 1, numberOfIntervals = 1000):
        cmd = "simulate(%s, stopTime=%f, numberOfIntervals=%d)"%(
            model,stopTime,numberOfIntervals)
        print self.__send__(cmd, timeout = 60*60*24)
    

    def simulate(self, model, stopTime=1, numberOfIntervals=1000):

        code="""
model simulate_$name
#set $sep = '' ##slurp
  extends ${name}(#slurp
#for key,value in $parameters.items() ##slurp
$sep$key = $value#slurp 
#set $sep = ", " ##slurp
#end for##slurp
);
end simulate_$name;
"""

        fd = NamedTemporaryFile(mode="w+", bufsize=0, prefix="simulate_")

        code = str(Template(code, searchList=[{
                            'name' : model,
                            'parameters' : self.parameters}]))   
        print code
        fd.write(code)

        self.loadFile(fd.name)

        fd.close()
        
        self.base_simulate("simulate_%s"%(model), stopTime, numberOfIntervals)

    def quit(self):
        self.shell.sendline("quit()")

    def list(self, name):
        print self.__send__("list(%s)"%(name))

    def val(self, name, time):
        return float(self.__send__("val(%s, %f)"%(name, time)).split("\n")[-2])

    def set_parameter(self, name, value):
        self.parameters[name] = value;

    def optimize_parameter(self, name, space, cost):
        if name not in self.opt_parameters:
            self.opt_parameters[name] = {
                'space': space,
                'cost': []}
            
        for c in cost:
            self.opt_parameters[name]["cost"].append({
                    'name': c[0],
                    'startCost': c[1],
                    'stopCost': c[2]
                    })


    def optimize_scipy(self, name, stopTime = 1, numberOfIntervals=1000, number=10, tolerance=0.01**2, step=0.01):

        # Define wrapper function to be used with SciPy
        def F(*args):
            print "Argumentos:\n\t", args

            obj = args[-4]
            name = args[-3]
            stopTime = args[-2]
            numberOfIntervals = args[-1]


            for i, param in zip(xrange(len(args[0])), obj.opt_parameters.keys()):
                obj.opt_parameters[param]['value'] = args[0][i]


            return obj.optimize_get_cost(
                name = name,
                stopTime = stopTime,
                numberOfIntervals = numberOfIntervals)


        # Search a "good" initial point.
        res = self.optimize(name, stopTime, numberOfIntervals, number)
        
        def callbk(xk):
            print "Xk: ", xk

        # return scipy.optimize.fmin_l_bfgs_b(
        #     F, 
        #     res[1][1:], 
        #     approx_grad = True,
        #     args = (self, name, stopTime, numberOfIntervals),
        #     iprint = 1)

        return scipy.optimize.fmin(
            F,
            res[1][1:],
            args = (self, name, stopTime, numberOfIntervals),
            retall = True,
            disp = True,
            callback = callbk)

    def optimize_get_cost(self, name, stopTime=1, numberOfIntervals=1000):

        code = """
model opt_$name
#set $sep = ''# #slurp
  extends ${name}(#slurp
#for key, value in $opt_parameters.items() ##slurp
$sep$key = $value['value']#slurp 
#set $sep = ", "##slurp
#end for#);

   
  Real opt_cost;
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
  Real opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']};
#end for##slurp
#end for#


equation
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
  der(opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']}) = if ((time < $v['startCost']) or (time > $v['stopCost'])) then 0 else $v['name'];
#end for##slurp
#end for#
  opt_cost = #slurp
#set $sep = ''##slurp
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
${sep}(opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']}/(${v['stopCost']}-${v['startCost']}))#slurp 
#set $sep = " + "##slurp
#end for##slurp
#end for#;
end opt_$name;
"""
            
        fd = NamedTemporaryFile(mode="w+", bufsize=0, prefix="Modelica_Optimize")
            
        code = str(Template(code, searchList=[{
                        'opt_parameters': self.opt_parameters,
                        'name' : name}]))
        
        print code
        
        fd.write(code)

        self.loadFile(fd.name)
        
        fd.close()
        
        self.base_simulate("opt_%s"%(name), stopTime, numberOfIntervals)
        
        print "Returning: ", self.val("opt_cost", stopTime)
        
        return self.val("opt_cost", stopTime)


    def optimize_simple(self, name, stopTime = 1, numberOfIntervals=1000, number=10, tolerance=0.01**2, step=0.01):
        # Search a "good" initial point.
        res = self.optimize(name, stopTime, numberOfIntervals, number)

        # Result of optimization
        result = res[0:2]
        
        n = len(result[0])-1
        
        g = [(2*random()-1) * x for x in [step]*n]

        # First parameters:
        for i in xrange(n):
            self.opt_parameters[result[0][i+1]]['value'] = result[1][i+1]

        # Update parameters
        last_opt = self.opt_parameters
        for i in xrange(n):
            self.opt_parameters[result[0][i+1]]['value'] *= 1 + g[i]

        error = res[1][0]

        count = 0

        while error > tolerance:
                
            code = """
model opt_$name
#set $sep = ''# #slurp
  extends ${name}(#slurp
#for key, value in $opt_parameters.items() ##slurp
$sep$key = $value['value']#slurp 
#set $sep = ", "##slurp
#end for#);

   
  Real opt_cost;
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
  Real opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']};
#end for##slurp
#end for#


equation
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
  der(opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']}) = if ((time < $v['startCost']) or (time > $v['stopCost'])) then 0 else $v['name'];
#end for##slurp
#end for#
  opt_cost = #slurp
#set $sep = ''##slurp
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
${sep}(opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']}/(${v['stopCost']}-${v['startCost']}))#slurp 
#set $sep = " + "##slurp
#end for##slurp
#end for#;
end opt_$name;
"""

            fd = NamedTemporaryFile(mode="w+", bufsize=0, prefix="Modelica_Optimize")
            
            code = str(Template(code, searchList=[{
                            'opt_parameters': self.opt_parameters,
                            'name' : name}]))

            print code
            fd.write(code)

            self.loadFile(fd.name)

            fd.close()

            self.base_simulate("opt_%s"%(name), stopTime, numberOfIntervals)

            opt = self.val("opt_cost", stopTime)


            if opt > result[-1][0]:
                # The changes is not in a descent direction.
                # Search a new direction, and try with it.
                
                g = [(2*random()-1) * x for x in [step]*n]

                print count, " x ", tuple([opt] + map(
                        lambda x:x['value'], 
                        map(
                            lambda x:self.opt_parameters[x],
                            result[0][1:])))
                
                self.opt_parameters = last_opt

            else:
                # Descent direction.
                result.append(
                    tuple([opt] + map(
                            lambda x: x['value'],
                            self.opt_parameters.values()
                            )
                          )
                    )

                print count, " - ", result[-1]
                count += 1

                print "\n".join(map(lambda x: " , ".join(map(str,x)),result))



            # Update parameters
            last_opt = self.opt_parameters
            for i in xrange(n):
                self.opt_parameters[result[0][i+1]]['value'] *= 1 + g[i]

                



                


    def optimize(self, name, stopTime = 1, numberOfIntervals=1000, number=10):
        # Define "number" random values in space and 
        # return "vble" value with best cost function.

        result = []

        for i in xrange(number):

            code = """
model opt_$name
#set $sep = ''# #slurp
  extends ${name}(#slurp
#for key, value in $opt_parameters.items() ##slurp
$sep$key = $value['value']#slurp 
#set $sep = ", "##slurp
#end for#);

   
  Real opt_cost;
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
  Real opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']};
#end for##slurp
#end for#


equation
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
  der(opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']}) = if ((time < $v['startCost']) or (time > $v['stopCost'])) then 0 else $v['name'];
#end for##slurp
#end for#
  opt_cost = #slurp
#set $sep = ''##slurp
#for key,value in $opt_parameters.items() ##slurp
#for v in $value['cost'] ##slurp
${sep}(opt_cost_${key}_${v['name']}_${v['startCost']}_${v['stopCost']}/(${v['stopCost']}-${v['startCost']}))#slurp 
#set $sep = " + "##slurp
#end for##slurp
#end for#;
end opt_$name;
"""

            fd = NamedTemporaryFile(mode="w+", bufsize=0, prefix="Modelica_Optimize")
            
            for par in self.opt_parameters:
                self.opt_parameters[par]['value'] = random() * (
                    self.opt_parameters[par]['space'][1] - self.opt_parameters[par]['space'][0]
                    ) + self.opt_parameters[par]['space'][0]

            code = str(Template(code, searchList=[{
                            'opt_parameters': self.opt_parameters,
                            'name' : name}]))

            print code
            fd.write(code)

            self.loadFile(fd.name)

            fd.close()

            self.base_simulate("opt_%s"%(name), stopTime, numberOfIntervals)

            opt = self.val("opt_cost", stopTime)

            result.append(
                tuple([opt] + map(
                        lambda x: x['value'],
                        self.opt_parameters.values()
                        )
                      )
                )

            print i, " - ", result[-1]
            
        result.sort(cmp=lambda x,y: 
                    ((x[0] > y[0]) and 1) or ((x[0] < y[0]) and  -1) or 0)
        return [["cost"] + self.opt_parameters.keys()] + result

            

