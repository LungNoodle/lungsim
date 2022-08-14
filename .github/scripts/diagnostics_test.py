from aether import diagnostics

initial_state = diagnostics.get_diagnostics_on()
diagnostics.set_diagnostics_on(True)
final_state = diagnostics.get_diagnostics_on()
print("Initial state: {0}, Final state: {1}".format(initial_state, final_state))

if initial_state != 0 or final_state != 1:
  print("Test failed.")
  sys.exit(1)
else:
  print("Test succeeded.")
