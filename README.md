# WHD model

## To Build

`g++` is the default compiler.


```
make
```
For more details please refer to the `Makefile`
Or compile diretly by typing:
./make.sh

## To Run

```
./build/bin/malaria_model -c config.json -o outputs
```

### Mandatory Commandline Options

- `-c [FILE]`, user specified simulation configuration file. A simplified (not complete) example `config_whd.json` is given below.

```
{
    "simulation" : {
        "schema_file_name" : "schema.config.json",
        "total_steps" : 60,
        "total_population" : 100
    },

    "whd" : {
        "whd_level" : 0,
        "II50" : 2000000000,
        "pmf" : 11.0,
        "gamma": 0.75,
                     
        "transition_matrix" : {
            "Jpk" : {
                "CL": 3.2500,
                "V2": 5.3750
            }
        }
    }
}
```

- `-o [OUTPUT_DIRECTORY]`

### Optional Commandline Options

- `-i [OUTPUT_PREFIX]`, if given will be used as prefix for all outputs of this run, otherwise a timestamp will be used to identify outputs. 


## To Test
We use the [Catch2](https://github.com/catchorg/Catch2) framework.

1. Build tests:
```
make test
Or compile diretly by typing:
./make_4_test.sh
```

2. Run all tests:
```
malaria_model$ ./build/bin/test
```

3. List all tests:
```
malaria_model$ ./build/bin/test -l
```

4. Run a specific test:
```
malaria_model$ ./build/bin/test [whd:test]
```

```
malaria_model$ ./build/bin/test [whd_multiple:test]
```

## Dependencies
C++ Libraries and the version installed on the development machine:

### .hpp Dependencies
These are header-only libraries that are included as part of this project. No installation is required.

- [Eigen](https://eigen.tuxfamily.org), v3.3.9, for linear algebra / matrix operations. Included in include/third_party directory.
- [Catch2](https://github.com/catchorg/Catch2), v2.1.2, a C++ test framework. `./test/third_party/`


## License

The WHD model is licensed under the GNU GENERAL PUBLIC LICENSE. See LICENSE for more information.

## Author

Babak Mahdavi Ardestani

babak.m.ardestani@gmail.com
