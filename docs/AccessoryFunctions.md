# AccessoryFunctions Overview

AccessoryFunctions contains miscellaneous methods that don't necessarily fit anywhere else. They may be useful.


### Download File

```python
from accessoryFunctions import accessoryFunctions
accessoryFunctions.download_file(address, output_name, hour_start=18, hour_end=6, day_start=5, day_end=6, timeout=600)
```

This method will allow for file downloads between specified hours (default 6PM to 6AM) or on specified days (default Saturday and Sunday). The `address` argument must the URL for a single file, not a folder.
If the download is interrupted by the end of the download window, it will be continued when the download window reopens. 

IMPORTANT: `output_name` must be an absolute path, and for whatever reason curl doesn't like the _~_ character, so it can't be used in this command.

### Make Path

```python
from accessoryFunctions import accessoryFunctions
accessoryFunctions.make_path(inpath)
```

Simple method to make a directory if it doesn't already exist. Only argument is the directory you want to create.

### Print Time

```python
from accessoryFunctions import accessoryFunctions
accessoryFunctions.printtime(string, start)
```

Prints out elapsed time since `start`, followed by the text contained in `string`. Prints it out in bold, so it looks good! 

### Dependency Check

```python
from accessoryFunctions import accessoryFunctions
present = accessoryFunctions.dependency_check(dependency)
```

Not sure if a program required in your script is installed and on the $PATH? Use this method to check for you - it will return `True` if the program is present, and `False` if it isn't.
A good way to warn users of your scripts that not everything needed is installed. 

NOTE: Doesn't check the version of the program, so don't rely on it for that. 

### Find Paired Reads

```python
from accessoryFunctions import accessoryFunctions
pairs = accessoryFunctions.find_paired_reads(directory, forward_id='_R1', reverse_id='_R2')
```

This method will look into the directory specified with the `directory` argument and return a list of files that appear to be read pairs. Default forward and reverse identifiers are '_R1' and '_R2', respectively, but these can be changed with the `forward_id` and `reverse_id` keyword arguments.

The return will be a nested list with the following format: [[forward_1, reverse_1], [forward_2, reverse_2]], and so on.

