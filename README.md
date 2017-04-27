# Fluigi
The repository that will host the fluigi binaries for the class


## Requirements

Please install the latest Java Runtime >1.8

## Usage
To run the program, use:
```
java [-Xmx????m] -jar fluigi filename [-i parameterFile] [-o outputformat] [-v|h]
```
The flags are as follows:

```
-i / --input : Input parameter file
-o / --output: output format
-v / --hdl : Verilog input
-h / --help : Prints out usage
```
The Java flag [-Xmx????m] will ensure that ???? Megabytes of RAM is allocated for the Java virtual machine. This might be necessary in scenarios where the computer might by default allocate a smaller chunk of memory.
