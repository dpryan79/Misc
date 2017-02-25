from docutils import nodes
from docutils.parsers.rst import Directive  # noqa
from sphinx.util.nodes import nested_parse_with_titles


class FooDirective(Directive):
    def run(self):
        """
        Return a section with children. These are (in order):
         - A title
         - An option list
         - A title
         - A section
           - A title
           - An option list
        Ideally the titles should be linkified in the final output.
        """
        env = self.state.document.settings.env

        sec2 = nodes.section(ids=['Commands:'])
        sec2.document = self.state.document
        sec2 += nodes.title('Commands:', 'Commands:')
        subsec = nodes.section(ids=['BED-file'])
        sec2 += subsec
        subsec += nodes.title('BED-file', 'BED-file')

        # Make the bottom-most option list
        n1 = nodes.option_list_item('',
             nodes.option_group('', nodes.option_string(text='--sub-opt1')),
             nodes.description('', nodes.paragraph(text='I am the help for the --sub-opt1 option')))
        n2 = nodes.option_list_item('',
             nodes.option_group('', nodes.option_string(text='--sub-opt2')),
             nodes.description('', nodes.paragraph(text='I am the help for the --sub-opt2 option')))
        #sub_optList = nodes.option_list('', *[n1, n2])
        subsec += nodes.option_list('', *[n1, n2])

        # Title/section for the sub-command
        #subsec = nodes.section(ids=['BED-file'])
        #subsec += nodes.title('BED-file', 'BED-file', title_level=2)
        #subsec += sub_optList

        # Title for all sub-commands
        #title = nodes.title('', text='Commands:')
        #subsec = nodes.section('', *[title, subsec], ids=['Commands:'])

        # A title for the option list
        #optTitle = nodes.title('', text='Options:')
        sec = nodes.section(ids=['Options:'])
        sec = nodes.title('Options:', 'Options:')

        # A more general option list
        n1 = nodes.option_list_item('',
             nodes.option_group('', nodes.option_string(text='--opt1')),
             nodes.description('', nodes.paragraph(text='I am the help for the --opt1 option')))
        n2 = nodes.option_list_item('',
             nodes.option_group('', nodes.option_string(text='--opt2')),
             nodes.description('', nodes.paragraph(text='I am the help for the --opt2 option')))
        nol = nodes.option_list('', *[n1, n2])

        return [sec, nol, sec2]


def setup(app):
    app.add_directive('foo', FooDirective)
    return {'version': '0.1'}
