---
layout: post
title:  "Docker Basic"
date:   2021-07-30 07:10:10 +0800
usemathjax: true
categories: jekyll update
---

## docker basics

### installation 
- install on ubuntu

```bash
sudo apt install docker.io
```

- check whether docker client and server is properly installed

```bash
docker version
```

### image and container

- image(镜像)相当于虚拟机模板或类
- container（容器）相当于虚拟机或类的实例

- check locally available images

```bash
docker images
# docker image ls will also work
```

- check running containers

```bash
docker container ls
```

- show both running and stopped container

```bash
docker container ls -a
```

- remove docker images

```bash
docker rmi image_id
# docker image rm image_id will also work
```

- remove container instances

```bash
docker container rm container_id
```

- run a exutable with docker container from a image

```bash
# run a exutable with docker container from a image
# docker container run name.or.id.of.image path.to.exutable parameters
docker run ubuntu echo 'hello world'
# 'hello world'
# if you want docker remove the container after exited, add --rm
# docker container --rm run name.or.id.of.image path.to.exutable
docker run --rm ubuntu echo 'hello world'
# 'hello world'
# If the docker container has a prespecified entry point
# docker run --rm image.with.entry.point parameters
docker run --rm jinyf1998/cowsay.v0 goodbye cold world
# ____________________
#< goodbye cold world >
# --------------------
#        \   ^__^
#         \  (oo)\_______
#            (__)\       )\/\
#                ||----w |
#                ||     ||


``` 

```bash
docker run -it ubuntu /bin/bash
```

```bash
docker run --rm 
```


### docker volume / external storage


### Dockerfile


### dockerhub


### Networking in docker


### Reference

- [The Docker Book](https://dockerbook.com/)
