#!/usr/bin/env python
#coding:utf-8

######################################################
#
# gitx - git speed up tool with set proxy
# written by LiuLiqiu
#
######################################################

import sys
import os
import commands
import platform
import re
import ConfigParser
from subprocess import call

import requests

from . import __version__

__OS__ = platform.platform()
__VERSION__ = __version__
__GIT_CONFIG_PATH__ = os.environ['HOME'] + '/.gitconfig'
__GITX_CONFIG_PATH__ = os.environ['HOME'] + '/.gitxconfig'

#parse gitx config file
def parse_gitx_config(gitx_config_path):
    config = ConfigParser.SafeConfigParser()
    try:
        config.read(gitx_config_path)
        return config
    except ConfigParser.ParsingError as e:
        print 'parsing errors!'
        raise e

#init gitx config file git path content
def gitx_init_git_path(path):
    config = ConfigParser.SafeConfigParser()
    config.read(__GITX_CONFIG_PATH__)
    try:
        config.add_section('git')
    except ConfigParser.DuplicateSectionError as e:
        pass
    config.set('git', 'path', path)
    config.write(open(__GITX_CONFIG_PATH__, 'w'))

#init gitx config file proxy content
def gitx_init_proxy(url):
    config = ConfigParser.SafeConfigParser()
    config.read(__GITX_CONFIG_PATH__)
    if '@' in url:
        pat = re.compile(r'(.*):(.*)@(.*):(.*)')
        r = pat.match(url.strip('http://'))
        if 4 != len(r.groups()):
            print 'proxy url error,please check the url address'
            return
        username = r.group(1)
        appkey = r.group(2)
        address = r.group(3)
        port = r.group(4)
        try:
            config.add_section('account')
        except ConfigParser.DuplicateSectionError as e:
            pass
        config.set('account', 'username', username)
        config.set('account', 'appkey', appkey)
    else:
        pat = re.compile(r'(.*):(.*)')
        r = pat.match(url.strip('http://'))
        if 2 != len(r.groups()):
            print 'proxy url error,please check the url address'
            return
        address = r.group(1)
        port = r.group(2)
        try:
            config.remove_section('account')
        except ConfigParser.Error as e:
            pass

    try:
        config.add_section('proxy')
    except ConfigParser.DuplicateSectionError as e:
        pass
    config.set('proxy', 'address', address)
    config.set('proxy', 'port', port)
    config.write(open(__GITX_CONFIG_PATH__, 'w'))

# find out whether git is installed on system
def install_git():
    select = raw_input('Do you want to install git now(y/N):')
    if select.lower() in ('yes', 'y'):
        OS = __OS__.lower()
        try:
            if ('debian' in OS) or ('ubuntu' in OS):
                call('sudo apt-get install -y git'.split())
            elif ('centos' in OS) or ('redhat' in OS):
                call('sudo yum install -y git'.split())
            elif 'darwin' in OS:
                call('sudo brew install git'.split())
            else:
                print "Can't find proper package for your system,please install by yourself!"
                exit(1)
        except OSError as e:
            print 'install command error!'
            raise e
    else:
        exit(0)

def find_git(git_path):
    # check the setting git path
    if git_path:
        status, out = commands.getstatusoutput('ls ' + git_path)
        if status:
            print 'There seem to be not git on the set path,please install git first!'
            install_git()
    # check git default path
    else:
        status1, out = commands.getstatusoutput('ls /usr/bin/git')
        status2, out = commands.getstatusoutput('ls /usr/libexec/git-core/git')
        status3, out = commands.getstatusoutput('ls /usr/lib/git-core/git')
        # git is not installed on system
        if status1 and status2 and status3:
            print 'There seem to be not git on your system,please install git first!'
            install_git()
        # git is installed,but git exec file not in /usr/bin/ directory
        elif status1 and not (status2 and status3):
            # git exec file in /usr/libexec/git-core/ directory
            if not status2:
                cp_command = ('sudo cp /usr/libexec/git-core/git /usr/bin/git').split()
                try:
                    call(cp_command)
                except OSError as e:
                    print 'cp command execute error!'
                    raise e
            # git exec file in /usr/lib/git-core/ directory
            elif not status3:
                cp_command = ('sudo cp /usr/lib/git-core/git /usr/bin/git')
                try:
                    call(cp_command)
                except OSError as e:
                    print 'cp command execute error!'
                    raise e

# execute git operation
def git_operation(argv, git_exec_path):
    str = git_exec_path + ' '
    for i in range(1, len(argv)):
        str = str + '\"' + argv[i] + '\" '
    try:
        os.system(str)
    except Exception,e:
        raise e

#get the username and appkey from config file
def get_username_appkey(config):
    username = None
    appkey = None
    if 'account' in config.sections():
        if 'username' in config.options('account'):
            username = config.get('account','username')
        if 'appkey' in config.options('account'):
            appkey = config.get('account','appkey')
        if username and appkey:
            return username + ':'+appkey
        else:
            return None
    else:
        return None

#get the real proxy server address
def get_proxy_server_address(config):
    if 'proxy' in config.sections():
        if 'address' in config.options('proxy') and 'port' in config.options('proxy'):
            real_proxy = config.get('proxy', 'address') + ':' + config.get('proxy', 'port')
        else:
            print 'no proxy server set in your config file!'
            real_proxy = None
        return real_proxy
    else:
        return None

#validate username and appkey
def validate_username_appkey(username_appkey,proxy_server):
    if proxy_server:
        validate_url = 'https://github.com'
        if username_appkey:
            proxy = 'http://'+username_appkey+'@'+proxy_server
        else:
            proxy = 'http://'+proxy_server
        proxies = { 'http':proxy, 'https': proxy, }
        try:
            result = requests.head(validate_url,proxies=proxies,timeout=30)
            if 200 == result.status_code:
                return proxy
        except requests.exceptions.RequestException as e:
            print 'validate failed or some other problem occured,not using proxy,more error information:'
            print e
    return 0

# set global proxy for git
def set_proxy(proxy_address,git_path):
    try:
        set_proxy_command = (git_path + ' config --global http.proxy '+proxy_address).split()
        call(set_proxy_command)
        print 'Set proxy successed!'
    except OSError as e:
        print 'Set proxy failed,more error information:'
        raise e
    try:
        instead_http_command = (git_path + ' config --global url.http://github.com/.insteadOf git://github.com/').split()
        call(instead_http_command)
        instead_https_command = (git_path + ' config --global url.https://github.com/.insteadOf git@github.com:').split()
        call(instead_https_command)
    except OSError as e:
        print 'stead git protocol failed,more error information:'
        raise e

# unset the proxy setting
def unset_proxy(git_path):
    unset_proxy_command = (git_path + ' config --global --unset http.proxy').split()
    call(unset_proxy_command)

# clone operation
def clone(argv,config,git_exec_path):
    address_list = argv[2:]
    address = None
    flag_tuple = (r'https://github.com', r'git@github.com', r'git://github.com')

    for temp in address_list:
        if not temp.startswith(r'-'):
            address = temp
            break

    if address == None:
        git_operation(argv, git_exec_path)

    elif address.startswith(flag_tuple) or (address.startswith(r"https://") and (r"@bitbucket.org" in address)):
        proxy_address = validate_username_appkey(get_username_appkey(config), get_proxy_server_address(config))
        if proxy_address:
            set_proxy(proxy_address, git_exec_path)
        git_operation(argv, git_exec_path)
        if proxy_address:
            unset_proxy(git_exec_path)

    elif address.startswith(r'git@bitbucket.org:'):
        proxy_address = validate_username_appkey(get_username_appkey(config), get_proxy_server_address(config))
        if proxy_address:
            set_proxy(proxy_address, git_exec_path)
        username = raw_input('Please input your bitbucket username:')
        command = git_exec_path + r' clone https://' + username + '@bitbucket.org/' + address.split(':')[1]
        try:
            call(command.split())
        except OSError as e:
            raise e
        if proxy_address:
            unset_proxy(git_exec_path)
    else:
        git_operation(argv, git_exec_path)

def command_line_runner():
    argv = sys.argv
    if 'proxy_set' == argv[1]:
        if len(argv) == 3:
            gitx_init_proxy(argv[2])
            exit(0)
        else:
            print 'proxy set error!'
            exit(1)
    elif 'path_set' == argv[1]:
        if len(argv) == 3:
            gitx_init_git_path(argv[2])
            exit(0)
        else:
            print 'git path set error!'
            exit(1)

    config = parse_gitx_config(__GITX_CONFIG_PATH__)
    if 'git' in config.sections() and 'path' in config.options('git'):
        git_path = config.get('git', 'path')
        git_exec_path = git_path
    else:
        git_path = None
        git_exec_path = '/usr/bin/git'

    find_git(git_path)

    if len(argv) < 2:
        git_operation(argv,git_exec_path)

    elif 'pull' == argv[1] or 'push' == argv[1]:
        proxy_address = validate_username_appkey(get_username_appkey(config), get_proxy_server_address(config))
        if proxy_address:
            set_proxy(proxy_address,git_exec_path)

        git_operation(argv,git_exec_path)

        if proxy_address:
            unset_proxy(git_exec_path)

    elif len(argv) == 2:
        if argv[1] == '--version':
            status,out = commands.getstatusoutput(git_exec_path + ' --version')
            print 'gitx version '+__VERSION__+'(based on ' +out+')'
        else:
            git_operation(argv,git_exec_path)

    elif 'clone' == argv[1]:
        clone(argv,config,git_exec_path)

    else:
        git_operation(argv,git_exec_path)

if __name__ == '__main__':
    command_line_runner()