{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
    :members:
    :inherited-members:

    {% if attributes %}
    .. rubric:: Attributes
    .. autosummary::
        {% for item in attributes %}            
        {%- if not item.startswith('_') %}
        {{ item }}
        {%- endif %}
        {% endfor %}
    {% endif %}

    {% if methods %}
    .. rubric:: Methods
    .. autosummary::
        {% for item in methods %}
        {%- if not item.startswith('_') or item in ['__call__'] %}
        {{ item }}
        {%- endif %}
        {% endfor %}
    {% endif %}
