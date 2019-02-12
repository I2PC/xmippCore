# xmippCore


**>>> To install Xmipp, please visit [this](https://github.com/I2PC/xmipp#xmipp) <<<**

However, if you still want to install **only** the XmippCore software, follow:

```
mkdir <xmipp-bundle>  # if it doesn't exist yet
cd <xmipp-bundle>
wget https://raw.githubusercontent.com/I2PC/xmipp/devel/xmipp -O xmipp
chmod 755 xmipp
./xmipp --help
mkdir src             # if not there yet
cd src
git clone https://github.com/I2PC/xmippCore/
cd - 
./xmipp config                # Configure compilation variables
./xmipp check_config:         # Check that the configuration is correct
./xmipp compile N xmippCore:  # Compile only xmippCore
./xmipp install [dir]         # Install at dir (by default, ./build)
```
