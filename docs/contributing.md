

`xolotl` is far from feature complete, and your contributions are welcome.

### Reporting Bugs


* Is it a bug? Are you sure the bug persists after you run `transpile` and `compile` and `xolotl.cleanup`?
* Describe what the expected behavior is, and what the actual behavior was


### Requesting Features

* Describe what you want
* Describe why you want it
* List papers that describe this mechanism, or original research that describes the feature you want

### Adding New Conductances/Synapses/Controllers

* Look at existing conductances/synapses/controllers and use them as a guideline
* If you're making a new conductance, put them in ``c++/conductances/<first_author_name>YY`` where `YY` is the two digit year
* Make sure you add a reference to the paper you're getting the conductance details from in a comment at the top of the file. 
* Run `xolotl.testConductances` to make sure your conductance file compiles correctly. 
* Run `x.show` on your conductance file to inspect the activation curves to make sure that it matches what you want it to do
* Send us a pull request
