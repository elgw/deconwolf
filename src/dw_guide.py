#!/usr/bin/env python

# https://stackoverflow.com/questions/14200721/how-to-create-a-menu-and-submenus-in-python-curses
# https://docs.python.org/3/library/curses.html
import curses
from curses import panel
import os
import re
import sys

# "pick" from
# https://github.com/wong2/pick/blob/master/pick/__init__.py


def do_nothing(*args):
    pass

def defaultLambda(c):
    c = c.lower()
    lambdas = {'dapi': 461,
                'cy5': 664,  # alexa 647
                'a594': 617,
                'ir800': 794,  # or 814 if it is alexa790
                'a700':723,
                'a488': 519,
                'tmr': 562}
    lam = lambdas[c]
    if lam is None:
        lam = 600
        print(f"Warning: don't know what wavelength that {c} has")
        print(f" setting it to {lam}, please set it correctly.")
    return lam


def getImFiles():
    ireg = '^.*\\.tiff?$'
    images = [f for f in os.listdir('./')
              if (os.path.isfile(os.path.join('./', f))
                  and (re.match(ireg, f) is not None)
                      and (re.match('dw\_', f) is None))]

    return images


def getChannel(imfile):
    m = re.match(r"(.*)\_[0-9]*\.tiff?", imfile)
    if m:
        # print(f"{m.group(1)}")
        return m.group(1)
    else:
        return None


def getChannels(imfiles):
    # Remove trailing .tif / .tiff
    # Remove trailing _ABC
    channels = []
    for imfile in imfiles:
        chan = getChannel(imfile)
        if chan is not None:
            channels.append(chan)

    channels = list(set(channels))
    return channels


KEYS_ENTER = (curses.KEY_ENTER, ord('\n'), ord('\r'))
KEYS_UP = (curses.KEY_UP, ord('k'))
KEYS_DOWN = (curses.KEY_DOWN, ord('j'))
KEYS_SELECT = (curses.KEY_RIGHT, ord(' '))

class Picker(object):
    """The :class:`Picker <Picker>` object
    :param options: a list of options to choose from
    :param title: (optional) a title above options list
    :param multiselect: (optional) if true its possible to select multiple values by hitting SPACE, defaults to False
    :param indicator: (optional) custom the selection indicator
    :param default_index: (optional) set this if the default selected option is not the first one
    :param options_map_func: (optional) a mapping function to pass each option through before displaying
    """

    def __init__(self, options, title=None, indicator='*', default_index=0, multiselect=False, multi_select=False, min_selection_count=0, options_map_func=None):

        if len(options) == 0:
            raise ValueError('options should not be an empty list')

        self.options = options
        self.title = title
        self.indicator = indicator
        self.multiselect = multiselect or multi_select
        self.min_selection_count = min_selection_count
        self.options_map_func = options_map_func
        self.all_selected = []

        if default_index >= len(options):
            raise ValueError('default_index should be less than the length of options')

        if multiselect and min_selection_count > len(options):
            raise ValueError('min_selection_count is bigger than the available options, you will not be able to make any selection')

        if options_map_func is not None and not callable(options_map_func):
            raise ValueError('options_map_func must be a callable function')

        self.index = default_index
        self.custom_handlers = {}

    def register_custom_handler(self, key, func):
        self.custom_handlers[key] = func

    def move_up(self):
        self.index -= 1
        if self.index < 0:
            self.index = len(self.options) - 1

    def move_down(self):
        self.index += 1
        if self.index >= len(self.options):
            self.index = 0

    def mark_index(self):
        if self.multiselect:
            if self.index in self.all_selected:
                self.all_selected.remove(self.index)
            else:
                self.all_selected.append(self.index)

    def get_selected(self):
        """return the current selected option as a tuple: (option, index)
           or as a list of tuples (in case multiselect==True)
        """
        if self.multiselect:
            return_tuples = []
            for selected in self.all_selected:
                return_tuples.append((self.options[selected], selected))
            return return_tuples
        else:
            return self.options[self.index], self.index

    def get_title_lines(self):
        if self.title:
            return self.title.split('\n') + ['']
        return []

    def get_option_lines(self):
        lines = []
        for index, option in enumerate(self.options):
            # pass the option through the options map of one was passed in
            if self.options_map_func:
                option = self.options_map_func(option)

            if index == self.index:
                prefix = self.indicator
            else:
                prefix = len(self.indicator) * ' '

            if self.multiselect and index in self.all_selected:
                format = curses.color_pair(1)
                line = ('{0} {1}'.format(prefix, option), format)
            else:
                line = '{0} {1}'.format(prefix, option)
            lines.append(line)

        return lines

    def get_lines(self):
        title_lines = self.get_title_lines()
        option_lines = self.get_option_lines()
        lines = title_lines + option_lines
        current_line = self.index + len(title_lines) + 1
        return lines, current_line

    def draw(self):
        """draw the curses ui on the screen, handle scroll if needed"""
        self.screen.clear()

        x, y = 1, 1  # start point
        max_y, max_x = self.screen.getmaxyx()
        max_rows = max_y - y  # the max rows we can draw

        lines, current_line = self.get_lines()

        # calculate how many lines we should scroll, relative to the top
        scroll_top = getattr(self, 'scroll_top', 0)
        if current_line <= scroll_top:
            scroll_top = 0
        elif current_line - scroll_top > max_rows:
            scroll_top = current_line - max_rows
        self.scroll_top = scroll_top

        lines_to_draw = lines[scroll_top:scroll_top+max_rows]

        for line in lines_to_draw:
            if type(line) is tuple:
                self.screen.addnstr(y, x, line[0], max_x-2, line[1])
            else:
                self.screen.addnstr(y, x, line, max_x-2)
            y += 1

        self.screen.refresh()

    def run_loop(self):
        while True:
            self.draw()
            c = self.screen.getch()
            if c in KEYS_UP:
                self.move_up()
            elif c in KEYS_DOWN:
                self.move_down()
            elif c in KEYS_ENTER:
                if self.multiselect and len(self.all_selected) < self.min_selection_count:
                    continue
                return self.get_selected()
            elif c in KEYS_SELECT and self.multiselect:
                self.mark_index()
            elif c in self.custom_handlers:
                ret = self.custom_handlers[c](self)
                if ret:
                    return ret

    def config_curses(self):
        try:
            # use the default colors of the terminal
            curses.use_default_colors()
            # hide the cursor
            curses.curs_set(0)
            # add some color for multi_select
            # @todo make colors configurable
            curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_WHITE)
        except:
            # Curses failed to initialize color support, eg. when TERM=vt100
            curses.initscr()

    def _start(self, screen):
        self.screen = screen
        self.config_curses()
        return self.run_loop()

    def start(self):
        return curses.wrapper(self._start)


def pick(*args, **kwargs):
    """Construct and start a :class:`Picker <Picker>`.
    Usage::
      >>> from pick import pick
      >>> title = 'Please choose an option: '
      >>> options = ['option1', 'option2', 'option3']
      >>> option, index = pick(options, title)
    """
    picker = Picker(*args, **kwargs)
    return picker.start()


class Menu(object):
    def __init__(self, config, stdscreen):
        self.window = stdscreen.subwin(0, 0)
        self.window.keypad(1)
        self.panel = panel.new_panel(self.window)
        self.panel.hide()
        self.config = config
        panel.update_panels()

        self.position = 0
        self.items = []
        self.updateItems()

    def updateItems(self):
        self.items = []

        self.items.append(("-- Settings (select with arrows and press <enter> to change)", 0, do_nothing))

        for name, value in config.items():
            self.items.append((f"{name} = {value}", 1, name))

        self.items.append(("-- Actions", 0, do_nothing))
        self.items.append(("run now!", "run"))
        self.items.append(("exit (and run later)", "exit"))
        self.items.append(("cancel", "cancel"))

    def navigate(self, n):
        self.position += n
        if self.position < 0:
            self.position = 0
        elif self.position >= len(self.items):
            self.position = len(self.items) - 1

    def display(self):
        self.panel.top()
        self.panel.show()
        self.window.clear()

        while True:
            self.window.refresh()
            curses.doupdate()
            for index, item in enumerate(self.items):
                if index == self.position:
                    mode = curses.A_REVERSE
                else:
                    mode = curses.A_NORMAL

                # msg = "%d. %s" % (index, item[0])
                msg = f"{item[0]}"
                self.window.addstr(1 + index, 1, msg, mode)

            key = self.window.getch()

            if key in [curses.KEY_ENTER, ord("\n")]:

                if self.position == len(self.items) - 1:
                    config['write'] = 0
                    config['run'] = 0
                    break

                if self.position == len(self.items) - 2:
                    config['write'] = 1
                    config['run'] = 0
                    break

                if self.position == len(self.items) - 3:
                    config['write'] = 1
                    config['run'] = 1
                    break;
                else:
                    item = self.items[self.position]
                    if(item[1] == 0):
                        # if type 0: run the function
                        item[2]()
                    if(item[1] == 1):
                        curses.echo()
                        # Move cursor to right of current
                        # read until a KEY_ENTER
                        posy = self.position+1
                        posx = len(item[0])+5
                        self.window.addstr(self.position+1, len(item[0])+2, "-> ", curses.A_NORMAL)
                        while not self.window.getch() in [curses.KEY_ENTER, ord("\n")]:
                            # just wait
                            0
                        cstring = self.window.instr(posy, posx, 10);
                        try:
                            value = str(float(cstring))
                            self.config[item[2]] = value;
                        except ValueError:
                            0

                        self.updateItems()
                        self.window.clear()

            elif key == curses.KEY_UP:
                self.navigate(-1)

            elif key == curses.KEY_DOWN:
                self.navigate(1)

        self.window.clear()
        self.panel.hide()
        panel.update_panels()
        curses.doupdate()


class MyApp(object):
    def __init__(self, stdscreen, config):
        self.screen = stdscreen
        curses.curs_set(0)

        main_menu = Menu(config, self.screen)
        main_menu.display()
        self.config = main_menu.config


def writeConfig():
    config['threads'] = int(round(float(config['threads'])))
    config['tilesize'] = int(round(float(config['tilesize'])))

    # Write out commands to generate PSF

    # Write out command to deonwolf
    print("Writing things to do to `dw_jobs`")
    with open('dw_jobs', "w") as file:
        for c in channels:
            lambd = config[c + '_nm']
            file.write(f"dw_bw --NA {config['NA']} "
                  f"--lambda {lambd} "
                  f"--ni {config['ni']} "
                  f"--threads {round(config['threads'])} "
                  f"--resxy {config['xy_res_nm']} "
                  f"--resz {config['z_res_nm']} "
                  f"PSF_{c}.tif "
                  f"\n")
        for f in imfiles:
            c = getChannel(f)
            if c is None:
                continue
            it = int(round(float((config[c + '_iter']))))
            psf = f"PSF_{c}.tif"
            file.write(f"dw --iter {it} "
                       f"--threads {config['threads']} "
                       f"--tilesize {config['tilesize']} "
                       f"{f} {psf}\n")


if __name__ == "__main__":
    #

    # Identify images
    imfiles = getImFiles()
    if len(imfiles) == 0:
        print(f"No files ending with .tif, or .tiff")
        print(f"Please run dw_guide in a folder with images")
        sys.exit(1)

    print(f"Found {len(imfiles)} images to be deconvolved")

    # Get defaults for known channels
    channels = getChannels(imfiles)

    config_BS1_100 = {'NA' : 1.45,
              'ni' : 1.15,
              'xy_res_nm' : 130,
              'z_res_nm' : 200
              }

    config_BS1_60 = {'NA' : 1.45,
              'ni' : 1.15,
              'xy_res_nm' : 130,
              'z_res_nm' : 200
              }

    config_BS2_40 = {'NA' : 1.4,
                     'ni' : 1.15,
                     'zy_res_nm' : 100,
                     'z_res_nm' : 600}


    configs = [config_BS1_100, config_BS1_60, config_BS2_40]
    option, index = pick(['BS1 100x', 'BS1 60x', 'BS2 40'], 'Please choose a profile')
    config = configs[index]

    config['threads'] = 10
    config['tilesize'] = 2048

    for c in channels:
        config[c + "_iter"] = 50
        config[c + "_nm"] = defaultLambda(c)

    print(f"Found {len(channels)} different channels/emitters")

    # Add config for discovered channels
    m = curses.wrapper(MyApp, config)
    if config['write']:
        writeConfig()

    if config['run']:
        print("running 'bash dw_jobs'")
        print("restart it any time")
        os.system('bash dw_jobs')
