"""Some pre-defined networks."""
from models.networks._base import Network, kon_sigmoid
from models.networks.toggle_switch import network as toggle_switch
from models.networks.repressilator import network as repressilator

__all__ = ['Network', 'kon_sigmoid', 'repressilator', 'toggle_switch']
