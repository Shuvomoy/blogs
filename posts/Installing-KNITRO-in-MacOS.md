@def title = "Installing KNITRO in MacOS"
@def published ="June 6, 2022"
@def tags =["programming", "Julia"]

# Installing KNITRO in MacOS

**Shuvomoy Das Gupta**

*June 6, 2022*

In this blog, we will discuss how to install KNITRO in MacOS. First, unzip the KNITRO.zip file by following the instructions [here](https://www.artelys.com/docs/knitro/1_introduction/installation/unix.html). Then copy the unzipped folder to your home directory. Location of the home directory can be found by pressing Command+Shift+H in finder. 

After that, run the following commands in terminal (.zsh).

```julia
nano ~/.zshrc 
```

This will open the `.zshrc` file in terminal. Add the following two lines in the file. (Keep in mind that the precise location will vary from user to user.)

```julia
export PATH="/Users/shuvomoy_das_gupta/knitro-13.0.1-MacOS-64/knitroampl:$PATH"

export DYLD_LIBRARY_PATH="/Users/shuvomoy_das_gupta/knitro-13.0.1-MacOS-64/lib:$DYLD_LIBRARY_PATH$"
```

Now press `Ctrl+O` to save and the press `Enter` to confirm the save. Then, quit `nano` by pressing `Ctrl+X`.

Now, source the `.zshrc` file by entering

```julia
source ~/.zshrc 
```

 in terminal.

Check if the path and library path has been loaded properly by running the following in the terminal.

```julia
echo $PATH      
```

```julia
echo $DYLD_LIBRARY_PATH
```

Copy the license file to the home folder.

While installing `KNITRO.jl` in `Julia`, if you have any library problem, then add the following lines in the `deps.jl` file located in `KNITRO/BuildNumber/deps/`.

```julia
const libknitro = "/Users/shuvomoy_das_gupta/knitro-13.0.1-MacOS-64/lib/libknitro.dylib"
const amplexe = "/Users/shuvomoy_das_gupta/knitro-13.0.1-MacOS-64/knitroampl"
```

