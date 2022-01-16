---
layout: post
title:  "Notes On Linux System Administration"
date:   2021-11-11 20:54:12 +0800
usemathjax: true
categories: jekyll update
---


### Check avaiable resources

- Avaliable memory: in `/proc/meminfo`

```bash
 free -h
#              total        used        free      shared  buff/cache   available
#Mem:            30G        5.3G         18G        120M        7.4G         25G
#Swap:          8.0G        2.0G        6.0G

```

### Daemon process

- Configuration files of systemd
  - Not modify
    - `/usr/lib/systemd/system`
    - `/lib/systemd/system`
  - Can be customized
    - `/etc/systemd/system`
    - `/run/systemd/system`

- List avaiable service

```bash
systemctl list-units --type=service
systemctl list-unit-files --type=service
```

- Show logging

```bash
journalctl
```

- <https://askubuntu.com/questions/903354/difference-between-systemctl-and-service-commands>


### ssh login without typing password

- Append local public key `~/.ssh/id_rsa.pub` to `~/.ssh/authorized_keys` in remote machine should work


### crontab

- Configuration: `minute(s) hour(s) day(s) month(s) weekday(s) command(s)`
- `/var/spool/cron/crontabs`
- See <https://kb.iu.edu/d/afiz>



### Account management
- Add account
  - `adduser`: a system command, not create home directory 
  - `useradd`: a perl wrapper for adduser, create home directory
- Remove account

```bash
userdel -r {user.name}
```

- Change user group

```bash
 usermod -a -G {group.name} {user.name}
```

### Network file system (NFS)


### NIS, LDAP and single sign on

- LDAP (Lightweight Directory Access Protocol)
- PAM (Pluggable Authentication Module system)

- <https://blog.csdn.net/developerinit/article/details/76141065>
- <https://www.cnblogs.com/lfdblog/p/9803276.html>

- Show all user with access

```bash
getent passwd
```

- Show password specific to LDAP

```bash
getent passwd --service=ldap
```


### Security related topics

- `/var/log/secure` provides useful loggings

- **Do not log in as root**

- <https://www.liquidweb.com/kb/how-do-i-set-up-setuid-setgid-and-sticky-bits-on-linux/>

