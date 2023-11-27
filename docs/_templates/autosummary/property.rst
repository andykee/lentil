:orphan:

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoattribute:: {{ fullname | replace("lentil.", "lentil::") }}

{# Normally this line would read
    .. auto{{ objtype }}:: {{ fullname | replace("lentil.", "lentil::") }}
but we've explicitly called autoattribute so that properties are rendered
in a way that is indistinguishable from attributes #}