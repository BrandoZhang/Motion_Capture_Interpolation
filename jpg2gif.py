import os
from PIL import Image

dir = input("Directory for jpeg series: ")
gif_name = input("GIF name to save: ")

duration = 17  # in ms
num_frames = 3261
frames = [(os.path.join(dir, f'pic{i:05}.jpeg'), duration) for i in range(num_frames)]
# Creates a list of PIL Image objects for each frame
images = [Image.open(frame[0]) for frame in frames]

# Saves the frames as an animated GIF
images[0].save(gif_name, save_all=True, append_images=images[1:], duration=[frame[1] for frame in frames], loop=0)

