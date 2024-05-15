import numpy as np
import json
import os

electrodes = np.zeros(32)

for i in range(4, 14):
    electrodes[i] = int(14 - i)

for i in range(20, 30):
    electrodes[i] = int(i - 9)

dac_channels = np.arange(1, 33)
dac_channels = np.array([str(channel) for channel in dac_channels])
mapping_dict = {}
print(len(dac_channels), len(electrodes))
for i, channel in enumerate(dac_channels):
    print(i)
    mapping_dict[channel] = electrodes[i]

print(mapping_dict)

name = 'default_artiq_config'
with open(str(name) + ".json", "w") as outfile: 
    json.dump(mapping_dict, outfile)
