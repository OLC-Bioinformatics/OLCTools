"""
CustomBox extends Box to provide:
- AttributeError for non-existent attributes
- Automatic creation of nested attributes
- Custom serialization to exclude unwanted keys and large values
- Additional functionalities: merge, key existence check, custom JSON
    serialization

Usage Example:
    >>> sample = CustomBox()
    >>> sample.mlst.profile = 23
    >>> sample.large_list = list(range(1000))
    >>> sample_dict = sample.to_dict(unwanted_keys=['large_list'])
    >>> print(sample_dict)
    {'mlst': {'profile': 23}}

    >>> other = {'virulence': {'profile': 45}}
    >>> sample.merge(other)
    >>> print(sample.virulence.profile)  # Output: 45

    >>> print(sample.key_exists('profile'))  # Output: True
    >>> print(sample.key_exists('non_existent'))  # Output: False

    >>> json_str = sample.to_json(unwanted_keys=['large_list'])
    >>> print(json_str)
    {
        "mlst": {
            "profile": 23
        },
        "virulence": {
            "profile": 45
        }
    }
"""

# Standard library imports
import json
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Union,
)

# Third-party imports
from box import Box


class CustomBox(Box):
    """
    CustomBox extends Box to provide:
    - AttributeError for non-existent attributes
    - Automatic creation of nested attributes
    - Custom serialization to exclude unwanted keys and large values
    - Additional functionalities: merge, key existence check, custom JSON
        serialization

    Usage Example:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> sample.large_list = list(range(1000))
        >>> sample_dict = sample.to_dict(unwanted_keys=['large_list'])
        >>> print(sample_dict)
        {'mlst': {'profile': 23}}

        >>> other = {'virulence': {'profile': 45}}
        >>> sample.merge(other)
        >>> print(sample.virulence.profile)  # Output: 45

        >>> print(sample.key_exists('profile'))  # Output: True
        >>> print(sample.key_exists('non_existent'))  # Output: False

        >>> json_str = sample.to_json(unwanted_keys=['large_list'])
        >>> print(json_str)
        {
            "mlst": {
                "profile": 23
            },
            "virulence": {
                "profile": 45
            }
        }
    """

    def __getattr__(self, item: str) -> Any:
        """
        Override __getattr__ to raise an AttributeError if the attribute
        does not exist.

        :param item: Attribute name
        :raises AttributeError: If the attribute does not exist
        :return: Attribute value

        Example usage:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> print(sample.mlst.profile)  # Output: 23
        >>> print(sample.mlst.non_existent)  # Raises AttributeError
        """
        if item not in self:
            raise AttributeError(f"No attribute '{item}'")
        return super().__getattr__(item)

    def set_attr(self, key: str, value: Any) -> None:
        """
        Custom method to set attributes, converting dictionaries to CustomBox
        objects.

        :param key: Attribute name
        :param value: Attribute value

        Example usage:
        >>> sample = CustomBox()
        >>> sample.set_attr('mlst', {'profile': 23})
        >>> print(sample.mlst.profile)  # Output: 23
        """
        if isinstance(value, dict) and not isinstance(value, Box):
            # Convert dictionaries to CustomBox objects
            value = CustomBox(value)
        # Directly set the attribute in the instance's __dict__
        self.__dict__[key] = value

    def to_dict(
        self,
        *,
        unwanted_keys: Optional[List[str]] = None,
        max_size: int = 1000
    ) -> Dict[str, Any]:
        """
        Convert the CustomBox object to a dictionary, excluding unwanted
        keys and large values.

        :param unwanted_keys: List of keys to exclude
        :param max_size: Maximum size for lists/dictionaries to include
        :return: Dictionary representation of the CustomBox object

        Example usage:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> sample.large_list = list(range(1000))
        >>> sample_dict = sample.to_dict(unwanted_keys=['large_list'])
        >>> print(sample_dict)
        {'mlst': {'profile': 23}}
        """
        unwanted_keys = unwanted_keys or []
        result = {}
        for key, value in self.items():
            # Skip unwanted keys
            if key in unwanted_keys:
                continue
            # Skip large lists/dictionaries
            if max_size and isinstance(value, (list, dict)) and \
               len(value) > max_size:
                continue
            # Recursively convert nested CustomBox objects
            if isinstance(value, Box):
                result[key] = value.to_dict(
                    unwanted_keys=unwanted_keys,
                    max_size=max_size
                )
            else:
                result[key] = value
        return result

    def merge(
        self,
        other: Union[Dict[str, Any], 'CustomBox']
    ) -> None:
        """
        Merge another dictionary or CustomBox into the current object.

        :param other: Dictionary or CustomBox to merge

        Example usage:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> other = {'virulence': {'profile': 45}}
        >>> sample.merge(other)
        >>> print(sample.virulence.profile)  # Output: 45
        """
        for key, value in other.items():
            if isinstance(value, dict):
                # If the key does not exist, create a new CustomBox
                if key not in self:
                    self[key] = CustomBox()
                # Recursively merge dictionaries
                self[key].merge(value)
            else:
                self[key] = value

    def key_exists(
        self,
        key: str
    ) -> bool:
        """
        Check if a key exists at any nested level.

        :param key: Key to check
        :return: True if the key exists, False otherwise

        Example usage:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> print(sample.key_exists('profile'))  # Output: True
        >>> print(sample.key_exists('non_existent'))  # Output: False
        """
        if key in self:
            return True
        # Recursively check nested CustomBox objects
        for value in self.values():
            if isinstance(value, CustomBox) and value.key_exists(key):
                return True
        return False

    def to_json(
        self,
        *,
        unwanted_keys: Optional[List[str]] = None,
        max_size: int = 1000
    ) -> str:
        """
        Serialize the CustomBox object to JSON, excluding unwanted keys
        and large values.

        :param unwanted_keys: List of keys to exclude
        :param max_size: Maximum size for lists/dictionaries to include
        :return: JSON string representation of the CustomBox object

        Example usage:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> sample.large_list = list(range(1000))
        >>> json_str = sample.to_json(unwanted_keys=['large_list'])
        >>> print(json_str)
        {
            "mlst": {
                "profile": 23
            }
        }
        """
        return json.dumps(
            self.to_dict(
                unwanted_keys=unwanted_keys,
                max_size=max_size
            ),
            indent=4
        )

    def to_file(
        self,
        file_path: str,
        *,
        unwanted_keys: Optional[List[str]] = None,
        max_size: int = 1000
    ) -> None:
        """
        Write the JSON dump of the CustomBox object to a file.

        :param file_path: Path to the file where the JSON will be written
        :param unwanted_keys: List of keys to exclude
        :param max_size: Maximum size for lists/dictionaries to include

        Example usage:
        >>> sample = CustomBox()
        >>> sample.mlst.profile = 23
        >>> sample.large_list = list(range(1000))
        >>> sample.to_file('metadata.json', unwanted_keys=['large_list'])
        """
        with open(file_path, 'w', encoding='utf-8') as outfile:
            json.dump(
                self.to_dict(
                    unwanted_keys=unwanted_keys,
                    max_size=max_size
                ),
                outfile,
                sort_keys=True,
                indent=4,
                separators=(',', ': ')
            )
