# PyManufSol
Framework for symbolically generating manufactured solutions satisfying prescribed conditions with applications to code verification.
This is the first version of the software that should be considered as a beta version as many improvements could be done. It is associated to the companion conference paper entitled "A GENERAL METHODOLOGY FOR SYMBOLICALLY GENERATING MANUFACTURED SOLUTIONS SATISFYING PRESCRIBED CONDITIONS - APPLICATION TO TWO-PHASE FLOWS EQUATIONS", presented at the VI International Conference on Numerical and Symbolic Computation: Developments and Applications (SYMCOMP 2023), March 30-31, 2023.  

# Usage
To run examples, go to the `Examples` directory
```bash
cd Examples
```
and run the following command:
```bash
python3 example.py
```
For example:
```bash
python3 SpinningBubbleExample.py
```
If an error is caught, you may need to add the following lines to the header of the python script:
```python
import os
import sympy as sp
import sys

script_path = os.path.dirname(os.path.realpath(__file__))
project_path = script_path+"/../."
sys.path.insert(0, project_path)
```
