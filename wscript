#!/usr/bin/env python3
# encoding: utf-8
# 2023 Harald Klimach <harald.klimach@dlr.de>

APPNAME = 'seeder'

top = '.'
out = 'build'

def options(opt):
    '''Building options provided by Musubi.
       Remember, all options can be displayed with waf --help.'''
    opt.recurse('bin')
    opt.recurse('aotus')
    opt.recurse('tem')
    opt.recurse('sdr')


def configure(conf):
    '''Project configuration'''
    import os
    conf.recurse('aotus', 'subconf')
    conf.recurse('bin', 'preconfigure')
    # Initialize the coco preprocessing tool
    conf.load('coco')
    conf.env['COCOSET'] = 'default.coco'
    if not conf.options.coco_reports:
      # Make coco silent, if not explicitly asked for reports:
      if conf.env.COCOFLAGS:
        conf.env.COCOFLAGS.insert(0, '-s')
        conf.env.COCOFLAGS.append('-ad')
      else:
        conf.env.COCOFLAGS = ['-s', '-ad']
    conf.recurse('tem')
    conf.recurse('polynomials')
    conf.recurse('sdr')
    conf.recurse('bin', 'postconfigure')


def build(bld):
    '''Build the Seeder project'''
    from revision_module import fill_revision_string
    bld.recurse('bin')
    if not (bld.cmd == 'docu' and bld.env.fordonline):
        bld.recurse('aotus')
    fill_revision_string(bld, subdir='sdr')
    bld(rule='cp ${SRC} ${TGT}', source=bld.env.COCOSET, target='coco.set')
    bld.recurse('tem')
    bld.recurse('polynomials')
    bld.recurse('sdr')

    if not (bld.cmd == 'docu'):
        bld.add_group()

        bld(
            rule = 'cp ${SRC} ${TGT[0].abspath()}',
            source = bld.path.find_or_declare('sdr/seeder'),
            target = 'seeder'
        )
        bld(
            rule = 'cp ${SRC} ${TGT[0].abspath()}',
            source = bld.path.find_or_declare('sdr/sdr_harvesting'),
            target = 'sdr_harvesting'
        )
        bld(
            rule = 'cp ${SRC} ${TGT[0].abspath()}',
            source = bld.path.find_or_declare('aotus/lua'),
            target = 'lua'
        )

#clean build directory and coco completely to create the build from scratch
def cleanall(ctx):
    from waflib import Options
    Options.commands = ['distclean'] + Options.commands
    ctx.exec_command('rm coco')
