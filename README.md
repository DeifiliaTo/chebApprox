# chebApprox

Implementation of chebyshev polynomial representation of functions


## Build and run

Install the [`stack`](https://docs.haskellstack.org/en/stable/README/) build tool, then:

```sh
stack build
```

and then one of:

```sh
stack run
stack exec [executable name] -- [command line arguments to executable]
```

To load the project into an interactive environment, one of:

```sh
stack ghci
stack ghci --no-load --ghci-options=[command line arguments to ghci]
```

## Install LLVM
### macOS

Using [Homebrew](https://brew.sh) on macOS:

```bash
brew install llvm-hs/llvm/llvm-8
```

### Debian/Ubuntu

For Debian/Ubuntu based Linux distributions, the [LLVM.org](https://llvm.org)
website provides binary distribution packages. Check
[apt.llvm.org](http://apt.llvm.org/) for instructions for adding the correct
package database for your OS version; for example:

```sh
echo "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial main" | sudo tee -a /etc/apt/sources.list
echo "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial-8 main" | sudo tee -a /etc/apt/sources.list
sudo apt-get update
```

And then:

```sh
apt-get install llvm-8-dev
```

