# Deployment
1. Build the image with ```cd Docker && docker build -t pintron .``` or download
   the most recent image from Docker Hub
2. Create the ```config``` directory, with two files: ```root_key.pub``` (the
ssh public key of the root) and  ```web_key.pub``` (the
ssh public key of the pintron user). The latter must contain all public keys of
the webservers that will be allowed to connect to this server.
3. Run container (the equivalent of starting a daemon) with
```docker run -d -P -v /home/user/PIntron-scripts/Docker/config:/config  --name "pintron-1" pintron```

# Design

The Dockerfile builds an image that contains a ssh server which is used to feed
the input data to pintron and to get back the results.
Therefore ports 80 and 22 must be open, or you must map the ssh port when
running the container.

PIntron is built inside the /home/app/pintron directory.
Input files are placed into the /home/app/input dir, while the results are
stored into /home/app/results.

If you want to map those directories to some specific locations on your server,
for instance because the input data are large, you can mount them as external
shares when running the image, with something like
```docker run -d -P --name "pintron-1" -v /media/HUGE_HARD_DISK:/home/app/input
pintron```. See
[Managing Data in Containers](http://docs.docker.com/userguide/dockervolumes/#volume-def)
for more information.

# Customization

If you want to track a version of PIntron different from stable (master), put
the desired version into ```environment/pintron_version```, when building the image.
