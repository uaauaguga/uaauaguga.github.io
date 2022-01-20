---
layout: post
title:  "Notes on v2ray Configuration"
date:   2022-01-16  10:42:56 +0800
usemathjax: true
categories: jekyll update
---

- Some notes on v2ray configuration, 以备不时之需


- Scripts for installation

```bash
https://git.io/v2ray.sh
```

- v2ray中的数据流向(摘自<https://toutyrater.github.io/>)。可以看出，直接交互的是客户端的outbound和服务器的inbound。

```bash
{浏览器} <--(socks)--> {V2Ray 客户端 inbound <-> V2Ray 客户端 outbound} <--(VMess)-->  {V2Ray 服务器 inbound <-> V2Ray 服务器 outbound} <--(Freedom)--> {目标网站}
```

- 几种出口协议
  - VMess
  - shadowsocks
  - freedom: 直连
  - blackhole: 只吞入数据,相当于阻止访问

### Explaination for the configuration file

- See <https://www.v2ray.com/en/configuration/overview.html>

- log, where to save logging information, here is an example, similar in both client and server

```json
 "log": {
    "loglevel": "warning",
    "access": "/var/log/v2ray/access.log",
    "error": "/var/log/v2ray/error.log"
  }
```


- For server side:

-  inbounds
- a list, each item is a dict, for each item, 
  - port: 提供v2ray服务的端口号
  - protocol: "vmess", 即v2ray的协议名称。另一个选项是"shadowsocks"
  - settings: a dict
    - clients: list客户端信息
    

- outbounds
  - a list


- routing
  - a dict


- For client side


- See <https://stackoverflow.com/questions/20778771/what-is-the-difference-between-0-0-0-0-127-0-0-1-and-localhost> difference between `0.0.0.0`, `localhost` and `127.0.0.1`


### 不同的client使用v2ray

#### 浏览器

#### 终端应用

- github

```bash
#这里的socks5://127.0.0.1:10808应当和本地v2ray配置文件的inbound一致
git config --global http.https://github.com.proxy socks5://127.0.0.1:10808
```

- curl

```bash
#这里的socks5h://localhost:10808应当和本地v2ray配置文件的inbound一致
curl -x socks5h://localhost:10808 www.google.com
```

- wget 默认不支持socks5代理...


### Tips

- 注意系统时间设置，如果VPS和本地机系统时间不一致v2ray无法正确使用

- Generate uuid in linux terminal : `uuidgen`


### Some useful resources

- <https://toutyrater.github.io/>

- <https://github.com/233boy/v2ray/wiki/V2Ray%E6%90%AD%E5%BB%BA%E8%AF%A6%E7%BB%86%E5%9B%BE%E6%96%87%E6%95%99%E7%A8%8B>

- <https://github.com/v2fly/v2ray-examples>