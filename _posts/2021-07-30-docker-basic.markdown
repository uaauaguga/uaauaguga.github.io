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

- Run docker without sudo

```bash
groupadd docker
usermod -aG docker $USER
newgrp docker
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

- start/stop/remove container instances

```bash
docker container start container_id
docker container stop container_id
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

- `docker container exec`和`docker container run`功能比较类似，区别在于


### docker volume (数据卷)
- external storage
- 有的时候你希望docker利用宿主机上的存储
- Docker volumes are directories that are not part of the container’s UFS

```bash
# 把当前路径下的programming目录mount到docker container 的/data路径下
# 如果想挂载多个目录,多次指定-v就可以了
docker run -it --rm -v $PWD/programming:/data ubuntu /bin/bash
```

### build customized docker image
- two ways
  - install dependency in docker interactive terminal, than commit a snapshot of the container with `docker commit`
  - using `Dockerfile`, and run `docker build`
- 优劣比较<https://zzq23.blog.csdn.net/article/details/80571262>
- The docker commit method is not currently recommended, as building with a Dockerfile is far more flexible and powerful
- 常用Dockerfile 关键字
  - `FROM`
  - `RUN`
  - `ENV`
  - `ADD`
  - `COPY`
  - `VOLUME`
  - `ENTRYPOINT`
  - `CMD`

- An example

```Dockerfile
FROM debian
MAINTAINER John Smith <john@smith.com>
RUN apt-get update && apt-get install -y cowsay fortune
COPY entrypoint.sh /
ENTRYPOINT ["/entrypoint.sh"]
```

- run

```bash
# Assume a Dockerfile is present in current directory
docker build -t target.image.name .
```


### docker hub
- pull docker image from docker hub with `docker pull`

- push your docker image to docker hub
  - tag the docker image to your own repo with `docker tag`
  - run `docker push`


### Reference

- [The Docker Book](https://dockerbook.com/)
- [Using Docker: Developing and Deploying Software with Containers](https://www.oreilly.com/library/view/using-docker/9781491915752/)
