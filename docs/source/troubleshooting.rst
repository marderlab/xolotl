.. highlight:: matlab


Troubleshooting
***************


On macOS, I get an annoying warning saying "xcrun: error: SDK "macosx10.13.4" cannot be located"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the following in your shell (not the MATLAB prompt)::

	sudo xcode-select -s /Applications/Xcode.app

On macOS, I get a warning saying that "Warning: Xcode is installed, but its license has not been accepted."
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, make sure you have XCode installed (not just the Command Line Tools -- the whole thing). You can get this from the Mac App Store. Then, open XCode and accept the license. You will have to do this only once.

I ran the quickstart, but I don't see anything
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Are you using a tiny screen? Some UI elements may go out of the frame on very small screens. To fix this, acquire the handle to the figure and change the position property. For example ::

  x.manipulate;
  manip = gcf;
  manip.Position = [100 100 34 56];

I get an error saying I don't have a compiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You need a C/C++ compiler. You need to follow MATLAB's instructions_ on how to get one, how to install one, and how to configure one. It may be helpful to also see our advice on compilers_.

.. _instructions: https://www.mathworks.com/support/compilers.html
.. _compilers: compilers.rst
