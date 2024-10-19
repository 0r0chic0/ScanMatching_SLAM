from ament_index_python import get_package_share_directory
from launch import LaunchDescription
from launch_ros.actions import Node
import os

def generate_launch_description():
    return LaunchDescription([
        Node(
            package='lab5_pkg',
            executable='scanmatch_node',
            name='scanmatch',
            # parameters=[config]
        )
    ])