These are some tests to try to figure out what is going on with linux builds.  In some cases (but not all), you may see an error like:

Traceback (most recent call last):
  File "simple_example.py", line 30, in <module>
    from PDSim.flow.flow import FlowPath
  File "CoolProp.pxd", line 11, in init PDSim.flow.flow
ValueError: CoolProp.CoolProp.State has the wrong size, try recompiling. Expected 88, got 136


and these docker environments are hoped to try to track down what is going on.
