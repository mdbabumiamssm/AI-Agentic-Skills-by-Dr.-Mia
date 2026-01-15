# React Best Practices Rules

## 1. File Structure & Naming

*   **PascalCase** for component filenames and directories: `components/UserProfile/UserProfile.tsx`.
*   **kebab-case** for utility files: `utils/date-format.ts`.
*   **Colocation**: Keep styles, tests, and types close to the component.
    ```text
    /UserProfile
      ├── UserProfile.tsx
      ├── UserProfile.test.tsx
      ├── UserProfile.module.css
      └── types.ts
    ```

## 2. Functional Components & Hooks

*   **Prefer Functional Components**: Avoid Class components unless absolutely necessary (e.g., Error Boundaries in older React versions).
*   **Typed Props**: Always define an interface for props.
    ```typescript
    interface ButtonProps {
      label: string;
      onClick: () => void;
      variant?: 'primary' | 'secondary';
    }
    ```
*   **Destructure Props**:
    ```typescript
    // BAD
    const Button = (props: ButtonProps) => <button>{props.label}</button>;

    // GOOD
    const Button = ({ label, onClick, variant = 'primary' }: ButtonProps) => ...
    ```

## 3. State Management

*   **Lift State Up**: Only when necessary.
*   **Avoid Derived State in State**: Compute values during render or use `useMemo`.
    ```typescript
    // BAD
    const [fullName, setFullName] = useState(firstName + ' ' + lastName);

    // GOOD
    const fullName = `${firstName} ${lastName}`;
    ```
*   **Use Context Sparingly**: For global themes or user sessions, not for high-frequency updates (to avoid re-renders).

## 4. Effects & Data Fetching

*   **No Data Fetching in useEffect**: Prefer libraries like **TanStack Query (React Query)** or **SWR**. They handle caching, loading states, and deduplication.
*   **Strict Dependency Arrays**: Never lie to the dependency array. If you need to omit something, likely your logic is wrong or you need `useEffectEvent` (experimental) or `useCallback`.

## 5. Performance

*   **Code Splitting**: Use `React.lazy` and `Suspense` for route-level components.
*   **Virtualization**: Use `react-window` or `virtuoso` for long lists.
*   **Stable References**: Wrap functions passed as props to optimized children in `useCallback`.

## 6. Accessibility (a11y)

*   **Semantic HTML**: Use `<button>`, `<header>`, `<nav>` instead of `<div>` with handlers.
*   **ARIA attributes**: Use only when semantic HTML is insufficient.
*   **Keyboard Navigation**: Ensure all interactive elements are focusable.

## 7. Testing

*   **React Testing Library**: Test behavior, not implementation details.
*   **User Events**: Use `userEvent` over `fireEvent` for more realistic interactions.
